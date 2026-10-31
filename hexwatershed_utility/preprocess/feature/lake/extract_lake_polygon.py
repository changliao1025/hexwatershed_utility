import logging
import os
from typing import Optional

from osgeo import ogr, osr
from pyearth.gis.gdal.gdal_vector_format_support import get_vector_driver_from_filename

logger = logging.getLogger(__name__)


def extract_lake_polygon(
    sFilename_in: str,
    sFilename_out: str,
    dThreshold_area: float = 5000.0,
    sFilter: Optional[str] = None,
    sField_area: Optional[str] = None,
    iFlag_exterior_only: bool = True,
) -> str:
    """
    Extract lake polygons exceeding a minimum area threshold from a vector dataset.

    Supports global (e.g. HydroLAKES, Natural Earth, GLWD) and regional
    (e.g. NHDPlus waterbodies) lake datasets in any GDAL-compatible vector format.

    Parameters
    ----------
    sFilename_in : str
        Path to input lake vector file (Shapefile, GeoJSON, GeoPackage, etc.).
    sFilename_out : str
        Path to output vector file (format determined by extension, e.g. GeoJSON).
    dThreshold_area : float, default=5000.0
        Minimum lake area threshold in km².
    sFilter : Optional[str], default=None
        SQL attribute filter expression (e.g. "LAKE = 1" or "Lake_type = 1").
    sField_area : Optional[str], default=None
        Optional field name containing pre-computed area in km² (e.g. "Lake_area", "AreaSqKm").
        If provided, uses this field for filtering before computing spherical area.
    iFlag_exterior_only : bool, default=True
        If True, extract only the exterior ring of the lake polygon (removing islands).
        If False, preserve interior rings.

    Returns
    -------
    str
        Path to the generated output vector file.
    """
    if not os.path.exists(sFilename_in):
        raise FileNotFoundError(f"Input lake file not found: {sFilename_in}")

    from pyearth.gis.geometry.calculate_polygon_area import calculate_polygon_area
    from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates

    pDataSource_in = ogr.Open(sFilename_in, 0)
    if pDataSource_in is None:
        raise RuntimeError(f"Could not open input vector dataset: {sFilename_in}")

    pLayer_in = pDataSource_in.GetLayer(0)
    if pLayer_in is None:
        raise RuntimeError(f"Could not get layer from dataset: {sFilename_in}")

    pSpatialRef_in = pLayer_in.GetSpatialRef()
    srs_target = osr.SpatialReference()
    srs_target.ImportFromEPSG(4326)

    # Determine coordinate transformation if needed
    oTransform = None
    if pSpatialRef_in is not None and not pSpatialRef_in.IsSame(srs_target):
        oTransform = osr.CoordinateTransformation(pSpatialRef_in, srs_target)

    # Apply attribute filter if provided
    if sFilter is not None:
        pLayer_in.SetAttributeFilter(sFilter)
    elif sField_area is not None:
        pLayer_in.SetAttributeFilter(f"{sField_area} >= {dThreshold_area}")

    # Prepare output dataset
    pDriver_out = get_vector_driver_from_filename(sFilename_out)
    if pDriver_out is None:
        raise ValueError(f"Could not determine vector driver for: {sFilename_out}")

    if os.path.exists(sFilename_out):
        pDriver_out.DeleteDataSource(sFilename_out)

    pDataSource_out = pDriver_out.CreateDataSource(sFilename_out)
    if pDataSource_out is None:
        raise RuntimeError(f"Could not create output datasource: {sFilename_out}")

    pLayer_out = pDataSource_out.CreateLayer("lake", srs_target, ogr.wkbPolygon)
    pLayer_out.CreateField(ogr.FieldDefn("lake_id", ogr.OFTInteger))
    pLayer_out.CreateField(ogr.FieldDefn("area", ogr.OFTReal))
    pLayerDefn_out = pLayer_out.GetLayerDefn()

    count = 0
    pLayer_in.ResetReading()
    for pFeature_in in pLayer_in:
        pGeometry_in = pFeature_in.GetGeometryRef()
        if pGeometry_in is None:
            continue

        pGeom = pGeometry_in.Clone()
        if oTransform is not None:
            pGeom.Transform(oTransform)

        eType = pGeom.GetGeometryType()
        # Handle both Polygon and MultiPolygon
        aPolygons = []
        if eType in (ogr.wkbPolygon, ogr.wkbPolygon25D):
            aPolygons.append(pGeom)
        elif eType in (ogr.wkbMultiPolygon, ogr.wkbMultiPolygon25D):
            for i in range(pGeom.GetGeometryCount()):
                aPolygons.append(pGeom.GetGeometryRef(i))

        for oPoly in aPolygons:
            if oPoly.IsEmpty():
                continue

            exterior_ring = oPoly.GetGeometryRef(0)
            if exterior_ring is None:
                continue

            if iFlag_exterior_only:
                pGeometry_out = ogr.Geometry(ogr.wkbPolygon)
                pGeometry_out.AddGeometry(exterior_ring.Clone())
            else:
                pGeometry_out = oPoly.Clone()

            aCoords_gcs = get_geometry_coordinates(pGeometry_out)
            if aCoords_gcs is None or len(aCoords_gcs) < 3:
                continue

            dArea = calculate_polygon_area(aCoords_gcs[:, 0], aCoords_gcs[:, 1]) / 1.0e6
            if dArea < dThreshold_area:
                continue

            pFeature_out = ogr.Feature(pLayerDefn_out)
            pFeature_out.SetField("lake_id", count + 1)
            pFeature_out.SetField("area", float(dArea))
            pFeature_out.SetGeometry(pGeometry_out)
            pLayer_out.CreateFeature(pFeature_out)
            count += 1

    pDataSource_out = None
    pDataSource_in = None

    logger.info("Extracted %d lake polygon(s) to %s", count, sFilename_out)
    return sFilename_out


if __name__ == "__main__":
    import sys

    if len(sys.argv) > 2:
        extract_lake_polygon(sys.argv[1], sys.argv[2])
