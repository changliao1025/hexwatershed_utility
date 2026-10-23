import json
import logging
import os
from pathlib import Path
from typing import Optional, Tuple

import netCDF4 as nc
import numpy as np
from pyearth.gis.geometry.calculate_polygon_area import calculate_polygon_area
from pyesm.mesh.unstructured.esm_create_unstructured_domain_file_full import (
    esm_create_unstructured_domain_file_full,
)

from hexwatershed_utility.postprocess.river.get_geometry import get_geometry

logger = logging.getLogger(__name__)


def convert_hexwatershed_json_to_river_parameter_file(
    sFilename_json_in: str,
    sFilename_river_parameter_in: str,
    sFilename_river_parameter_out: str,
    sFilename_river_domain_out: str,
    pWidth_in: Optional[float] = None,
    pDepth_in: Optional[float] = None,
) -> Tuple[str, str]:
    """
    Convert HexWatershed JSON outputs to river model (MOSART) parameter and domain NetCDF files.

    Parameters
    ----------
    sFilename_json_in : str
        Input HexWatershed JSON file.
    sFilename_river_parameter_in : str
        Template river parameter NetCDF file (for format compatibility).
    sFilename_river_parameter_out : str
        Target river parameter NetCDF output path.
    sFilename_river_domain_out : str
        Target river domain NetCDF output path.
    pWidth_in : float, optional
        River width scaling parameter.
    pDepth_in : float, optional
        River depth scaling parameter.

    Returns
    -------
    Tuple[str, str]
        (path_to_parameter_file, path_to_domain_file)
    """
    # Resolve relative paths to absolute paths
    if not os.path.isabs(sFilename_json_in):
        sFilename_json_in = os.path.abspath(sFilename_json_in)
    if not os.path.isabs(sFilename_river_parameter_in):
        sFilename_river_parameter_in = os.path.abspath(sFilename_river_parameter_in)
    if not os.path.isabs(sFilename_river_parameter_out):
        sFilename_river_parameter_out = os.path.abspath(sFilename_river_parameter_out)
    if not os.path.isabs(sFilename_river_domain_out):
        sFilename_river_domain_out = os.path.abspath(sFilename_river_domain_out)

    # open json and read data
    aID = list()
    aDnID = list()
    aCellID = list()
    aCellID_downslope = list()
    aElevation = list()
    aDrainage = list()
    aSlope = list()
    aFlowline_length = list()
    aLongitude = list()
    aLatitude = list()
    aArea = list()
    aDomainfrac = list()
    aFdir = list()
    aGxr = list()

    aHslp = list()
    aID = list()
    aLat = list()
    aLatixy = list()
    aLon = list()
    aLongxy = list()
    aNh = list()
    aNr = list()
    aNt = list()
    aRdep = list()
    aRlen = list()
    aRslp = list()
    aRwid = list()
    aRwid0 = list()
    aTslp = list()
    aTwid = list()
    aLatV_region = list()
    aLonV_region = list()
    with open(sFilename_json_in) as json_file:
        data = json.load(json_file)
        ncell = len(data)
        lID = 1
        for i in range(ncell):
            pcell = data[i]
            aID.append(lID)
            lID = lID + 1

            # aArea.append(dArea)
            aLongitude.append(float(pcell["dLongitude_center_degree"]))
            aLatitude.append(float(pcell["dLatitude_center_degree"]))
            aCellID.append(int(pcell["lCellID"]))
            aCellID_downslope.append(int(pcell["lCellID_downslope"]))
            aElevation.append(float(pcell["dElevation"]))
            dDrainage = float(pcell["dDrainage_area"])
            aDrainage.append(dDrainage)
            # qa for slope
            dSlope = float(pcell["dSlope_between"])
            if dSlope < 0.0005:
                dSlope = 0.0005
                pass

            aSlope.append(dSlope)
            dFlowline_length = float(pcell["dLength_flowline"])
            aFlowline_length.append(dFlowline_length)
            aDomainfrac.append(1.0)
            aFdir.append(1)

            aNh.append(0.1)
            aNr.append(0.05)
            aNt.append(0.05)
            aTwid.append(10.0)  # this parameter should be updated

            dummy_vertex = pcell["vVertex"]
            aVertex_lon = np.full(9, -9999, float)
            aVertex_lat = np.full(9, -9999, float)
            nVertex = len(dummy_vertex)
            for j in range(nVertex):
                aVertex_lon[j] = dummy_vertex[j]["dLongitude_degree"]
                aVertex_lat[j] = dummy_vertex[j]["dLatitude_degree"]

            aLonV_region.append(aVertex_lon)
            aLatV_region.append(aVertex_lat)
            dArea = float(pcell["dArea"])
            aLongitude_dummy = aVertex_lon[np.where(aVertex_lon != -9999)]
            aLatitude_dummy = aVertex_lat[np.where(aVertex_lat != -9999)]
            # dArea = calculate_polygon_area(aLongitude_dummy, aLatitude_dummy)
            # float(pcell['dArea'])
            aArea.append(dArea)
            # dGxr = dArea / dFlowline_length
            dGxr = dFlowline_length / dArea
            aGxr.append(dGxr)

        pass

    aLongitude_in = np.array(aLongitude).reshape(ncell)
    aLatitude_in = np.array(aLatitude).reshape(ncell)
    nCell = len(aCellID)
    aCellID = np.array(aCellID).reshape(ncell)
    aCellID_downslope = np.array(aCellID_downslope).reshape(ncell)
    aArea = np.array(aArea).reshape(ncell)

    aRwid, aRdep, aFlood_2yr_out = get_geometry(
        aLongitude_in,
        aLatitude_in,
        aCellID,
        aCellID_downslope,
        aArea,
        pWidth_in=pWidth_in,
        pDepth_in=pDepth_in,
    )

    # convert to numpy array

    aAreaTotal2 = np.array(aDrainage)
    aAreaTotal = aAreaTotal2

    aID = np.array(aID)

    aLonV_region = np.array(aLonV_region)  # .T
    aLatV_region = np.array(aLatV_region)  # .T

    # Build cell ID to index mapping for O(1) lookups instead of O(N) np.where scans
    cell_id_to_id = {int(cid): int(aid) for cid, aid in zip(aCellID, aID)}
    aDnID = [
        cell_id_to_id.get(int(down_id), -9999) if down_id != -1 else -9999
        for down_id in aCellID_downslope
    ]

    aDnID = np.array(aDnID)
    aDomainfrac = np.array(aDomainfrac)
    aElevation = np.array(aElevation)
    aRlen = np.array(aFlowline_length)
    aLon = np.array(aLongitude)
    aLat = np.array(aLatitude)
    aLongxy = aLon
    aLatixy = aLat

    aDomainfrac = np.array(aDomainfrac)
    aFdir = np.array(aFdir)
    aGxr = np.array(aGxr)
    aHslp = np.array(aSlope)
    aNh = np.array(aNh)
    aNr = np.array(aNr)
    aNt = np.array(aNt)
    aRdep = np.array(aRdep)
    aRwid = np.array(aRwid)
    aTwid = np.array(aTwid) * 5
    # fix tributary width

    aRslp = aHslp
    aTslp = aHslp

    # open mpas netcdf to get format
    with nc.Dataset(sFilename_river_parameter_in, "r") as pDatasets_in:
        netcdf_format = pDatasets_in.file_format

    aVariable = [
        "area",
        "areaTotal",
        "areaTotal2",
        "dnID",
        "domainfrac",
        "fdir",
        "frac",
        "gxr",
        "hslp",
        "CellID",
        "ID",
        "lat",
        "latixy",
        "lon",
        "longxy",
        "nh",
        "nr",
        "nt",
        "rdep",
        "rlen",
        "rslp",
        "rwid",
        "rwid0",
        "tslp",
        "twid",
    ]

    aUnit = [
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
    ]
    aLongName = [
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
    ]

    # write data out
    aData_out = [
        aArea,
        aAreaTotal,
        aAreaTotal2,
        aDnID,
        aDomainfrac,
        aFdir,
        aFdir,
        aGxr,
        aHslp,
        aCellID,
        aID,
        aLat,
        aLatixy,
        aLon,
        aLongxy,
        aNh,
        aNr,
        aNt,
        aRdep,
        aRlen,
        aRslp,
        aRwid,
        aRwid * 5.0,
        aTslp,
        aTwid,
    ]

    # now export
    nVariable = len(aVariable)
    aDimension_tuple = ("gridcell",)

    with nc.Dataset(sFilename_river_parameter_out, "w", format=netcdf_format) as pDatasets_out:
        pDatasets_out.createDimension("gridcell", nCell)
        for i in range(nVariable):
            pVar = pDatasets_out.createVariable(
                aVariable[i], float, aDimension_tuple, fill_value=-9999
            )
            aData_variable = aData_out[i]
            pVar[:] = aData_variable
            pVar.setncatts({"_FillValue": -9999})
            pVar.setncatts({"unit": aUnit[i]})
            pVar.setncatts({"long name": aLongName[i]})

    # create domain file

    aShape = aLon.shape
    if len(aShape) == 1:
        nrow = aShape[0]
        ncolumn = 1
    else:
        nrow = aShape[0]
        ncolumn = aShape[1]

    aLon.shape = (nrow, 1)
    aLat.shape = (nrow, 1)

    nrow, nvertex = aLonV_region.shape
    aLonV_region.shape = (nrow, 1, nvertex)
    aLatV_region.shape = (nrow, 1, nvertex)

    esm_create_unstructured_domain_file_full(
        aLon,
        aLat,
        aLonV_region,
        aLatV_region,
        sFilename_river_domain_out,
        aArea_in=aArea,
    )

    return sFilename_river_parameter_out, sFilename_river_domain_out
