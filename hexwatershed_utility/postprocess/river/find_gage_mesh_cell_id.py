import copy
import logging
import os
import sys

import netCDF4 as nc
import numpy as np
from osgeo import gdal, ogr, osr
from rtree.index import Index as RTreeindex

logger = logging.getLogger(__name__)


def find_gage_mesh_cell_id(
    aSitename_in,
    aLongitude_gage_in,
    aLatitude_gage_in,
    aDrainage_area_in,
    sFilename_domain_in,
    sFilename_parameter_in,
    iSearch_radius_in=None,
    dThreshold_drainage_in=None,
    dThreshold_difference_in=0.10,
    iFlag_data_km_in=0,
    dBuffer_in=0.3,
):
    """
    Snap river stream gauges to mesh cells using spatial indexing (R-tree)
    and drainage area matching.

    Parameters
    ----------
    aSitename_in : list or array
        Stream gauge site identifiers.
    aLongitude_gage_in : list or array
        Gauge longitudes (degrees).
    aLatitude_gage_in : list or array
        Gauge latitudes (degrees).
    aDrainage_area_in : list or array
        Gauge observed drainage areas (km2 or m2 depending on flag).
    sFilename_domain_in : str
        Path to river domain NetCDF file (containing xv, yv vertex coordinates).
    sFilename_parameter_in : str
        Path to river parameter NetCDF file (containing ID, dnID, areaTotal2).
    iSearch_radius_in : int, optional
        Optional search radius.
    dThreshold_drainage_in : float, optional
        Minimum drainage area threshold. Defaults to 100 km2.
    dThreshold_difference_in : float, default=0.10
        Maximum allowed fractional difference between observed and model drainage area.
    iFlag_data_km_in : int, default=0
        Flag indicating if parameter drainage area is already in km2 (1) or m2 (0).
    dBuffer_in : float, default=0.3
        Bounding box search buffer around gauge in degrees.

    Returns
    -------
    tuple of (list, list, list)
        (matched_gauge_indices, matched_cell_ids, matched_model_drainage_areas)
    """
    if dThreshold_drainage_in is not None:
        dThreshold_drainage = dThreshold_drainage_in / 1.0e6
    else:
        dThreshold_drainage = 1.0e2

    sFilename_geojson_mesh = sFilename_parameter_in.replace(".nc", "_mesh.geojson")
    if os.path.exists(sFilename_geojson_mesh):
        os.remove(sFilename_geojson_mesh)

    with nc.Dataset(sFilename_domain_in, "r") as pDatasets_domain:
        aXV = pDatasets_domain.variables["xv"][:].data
        aYV = pDatasets_domain.variables["yv"][:].data
        aXC = pDatasets_domain.variables["xc"][:].data

    with nc.Dataset(sFilename_parameter_in, "r") as aDatasets:
        for sKey, aValue in aDatasets.variables.items():
            if sKey == "CellID":
                aCellID = (aValue[:]).data
            if sKey == "ID":
                aID = (aValue[:]).data
            if sKey == "dnID":
                aDnID = (aValue[:]).data
            if sKey == "fdir":
                aFdir = (aValue[:]).data
            if sKey == "latixy":
                aLatitude = (aValue[:]).data
            if sKey == "longxy":
                aLongitude = (aValue[:]).data
            if sKey == "areaTotal2":
                if iFlag_data_km_in == 1:
                    aAccu = (aValue[:]).data
                else:
                    aAccu = (aValue[:]).data / 1.0e6  # from m2 to km2

    index_cell = RTreeindex()

    pDriver = ogr.GetDriverByName("GeoJSON")
    pDataset = pDriver.CreateDataSource(sFilename_geojson_mesh)
    pSpatial_reference_gcs = osr.SpatialReference()
    pSpatial_reference_gcs.ImportFromEPSG(4326)

    pLayer = pDataset.CreateLayer("cell", pSpatial_reference_gcs, ogr.wkbPolygon)
    pLayer.CreateField(ogr.FieldDefn("id", ogr.OFTInteger64))
    pLayer.CreateField(ogr.FieldDefn("cellid", ogr.OFTInteger64))
    pLayer.CreateField(ogr.FieldDefn("dnID", ogr.OFTInteger64))
    pLayer.CreateField(ogr.FieldDefn("drain", ogr.OFTReal))

    pLayerDefn = pLayer.GetLayerDefn()
    pFeature = ogr.Feature(pLayerDefn)
    nCell = len(aID)
    aBoundary = list()
    pGeometry_union = ogr.Geometry(ogr.wkbPolygon)

    for i in np.arange(nCell):
        ring = ogr.Geometry(ogr.wkbLinearRing)
        pXV0 = aXV[i]
        pXV = pXV0[pXV0 != -9999]
        nVertex = len(pXV)
        for j in range(nVertex):
            x1 = aXV[i, 0, j]
            y1 = aYV[i, 0, j]
            ring.AddPoint(x1, y1)

        ring.AddPoint(aXV[i, 0, 0], aYV[i, 0, 0])
        pPolygon = ogr.Geometry(ogr.wkbPolygon)
        pPolygon.AddGeometry(ring)
        aBoundary.append(pPolygon)
        pFeature.SetGeometry(pPolygon)
        lID = aID[i]

        pFeature.SetField("id", int(lID))
        pFeature.SetField("dnID", int(aDnID[i]))
        pFeature.SetField("drain", float(aAccu[i]))
        pLayer.CreateFeature(pFeature)

        left, right, bottom, top = pPolygon.GetEnvelope()
        pBound = (left, bottom, right, top)
        index_cell.insert(i, pBound)

        pGeometry_union = pGeometry_union.Union(pPolygon)

    pDataset = pLayer = pFeature = None

    nGage = len(aLongitude_gage_in)
    sFilename_geojson_gage = sFilename_parameter_in.replace(".nc", "_gage.geojson")
    if os.path.exists(sFilename_geojson_gage):
        os.remove(sFilename_geojson_gage)

    pDataset_gage = pDriver.CreateDataSource(sFilename_geojson_gage)
    pLayer_gage = pDataset_gage.CreateLayer(
        "gauge", pSpatial_reference_gcs, ogr.wkbPoint
    )
    pLayer_gage.CreateField(ogr.FieldDefn("name", ogr.OFTString))
    pLayer_gage.CreateField(ogr.FieldDefn("lon", ogr.OFTReal))
    pLayer_gage.CreateField(ogr.FieldDefn("lat", ogr.OFTReal))
    pLayer_gage.CreateField(ogr.FieldDefn("drain", ogr.OFTReal))
    pLayer_gage.CreateField(ogr.FieldDefn("drain1", ogr.OFTReal))
    pLayer_gage.CreateField(ogr.FieldDefn("qa", ogr.OFTReal))
    pLayer_gage.CreateField(ogr.FieldDefn("cellid", ogr.OFTInteger64))
    pLayerDefn = pLayer_gage.GetLayerDefn()

    aIndex_out = list()
    aCellID_out = list()
    aDrainage_area_out = list()

    for i in range(nGage):
        sName = aSitename_in[i]
        dLongitude_gage = float(aLongitude_gage_in[i])
        dLatitude_gage = float(aLatitude_gage_in[i])
        dDrainage_area = float(aDrainage_area_in[i])

        if np.isnan(dDrainage_area):
            continue

        pPoint = ogr.Geometry(ogr.wkbPoint)
        pPoint.AddPoint(dLongitude_gage, dLatitude_gage)
        if not pPoint.Within(pGeometry_union):
            continue

        dBuffer = dBuffer_in
        left = dLongitude_gage - dBuffer
        right = dLongitude_gage + dBuffer
        bottom = dLatitude_gage - dBuffer
        top = dLatitude_gage + dBuffer
        pBound = (left, bottom, right, top)

        aIntersect = list(index_cell.intersection(pBound))
        if len(aIntersect) == 0:
            continue

        nIntersect = len(aIntersect)
        dDiff_min = np.inf
        dIndex_min = -1
        dDrainage_mesh_cell_min = 0.0

        for j in range(nIntersect):
            lIndex = int(aIntersect[j])
            dDrainage_mesh_cell = float(aAccu[lIndex])
            lDownID = int(aDnID[lIndex])

            if lDownID < 0:
                continue

            dDiff_new = float(np.abs(dDrainage_mesh_cell - dDrainage_area))
            if dDiff_new < dDiff_min:
                dDiff_min = float(dDiff_new)
                dIndex_min = lIndex
                dDrainage_mesh_cell_min = dDrainage_mesh_cell

        if dIndex_min < 0:
            continue

        pBoundary = aBoundary[dIndex_min]
        lCellID = int(aID[dIndex_min])

        dRatio = float(
            dDiff_min / np.max([dDrainage_area, dDrainage_mesh_cell_min])
        )
        iFlag_qa = float(dRatio)

        if (
            dRatio < dThreshold_difference_in
            and dDrainage_area > dThreshold_drainage
        ):
            aIndex_out.append(i)
            aCellID_out.append(lCellID)
            aDrainage_area_out.append(dDrainage_mesh_cell_min)
        else:
            logger.debug(
                "Drainage area difference %s exceeds threshold %s for gauge %s",
                dRatio,
                dThreshold_difference_in,
                sName,
            )

        pFeature = ogr.Feature(pLayerDefn)
        pFeature.SetGeometry(pPoint)
        pFeature.SetField("name", str(sName))
        pFeature.SetField("lon", dLongitude_gage)
        pFeature.SetField("lat", dLatitude_gage)
        pFeature.SetField("drain", dDrainage_area)
        pFeature.SetField("drain1", dDrainage_mesh_cell_min)
        pFeature.SetField("qa", iFlag_qa)
        pFeature.SetField("cellid", lCellID)
        pLayer_gage.CreateFeature(pFeature)

    pDataset_gage = pLayer_gage = pFeature = None
    return aIndex_out, aCellID_out, aDrainage_area_out
