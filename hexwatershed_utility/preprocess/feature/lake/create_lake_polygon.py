import logging
import os
from osgeo import ogr, osr

logger = logging.getLogger(__name__)


def create_lake_polygon(
    sFilename_in: str,
    sFilename_out: str,
    dThreshold_area: float = 5000.0,
) -> str:
    """
    Extract large lake polygons from a shapefile (e.g. HydroLAKES) to GeoJSON.

    Parameters
    ----------
    sFilename_in : str
        Input lake shapefile path.
    sFilename_out : str
        Output GeoJSON path.
    dThreshold_area : float, default=5000.0
        Minimum lake area threshold in km2.

    Returns
    -------
    str
        Path to output GeoJSON file.
    """
    if not os.path.exists(sFilename_in):
        raise FileNotFoundError(f"Input lake file not found: {sFilename_in}")

    from pyearth.gis.geometry.calculate_polygon_area import calculate_polygon_area
    from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates

    pDriver_shp = ogr.GetDriverByName("ESRI Shapefile")
    pDriver_json = ogr.GetDriverByName("GeoJSON")

    pDataSource_in = pDriver_shp.Open(sFilename_in, 0)
    pLayer_in = pDataSource_in.GetLayer()
    pSpatialRef_in = pLayer_in.GetSpatialRef()
    srs_wkt = pSpatialRef_in.ExportToWkt()

    if os.path.exists(sFilename_out):
        pDriver_json.DeleteDataSource(sFilename_out)

    pDataSource_out = pDriver_json.CreateDataSource(sFilename_out)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    pLayer_out = pDataSource_out.CreateLayer("lake", srs, ogr.wkbPolygon)
    pField = ogr.FieldDefn("area", ogr.OFTReal)
    pLayer_out.CreateField(pField)
    pLayerDefn = pLayer_out.GetLayerDefn()
    pFeature_out = ogr.Feature(pLayerDefn)

    sProjection = srs.ExportToWkt()
    iProjection = 1 if srs_wkt != sProjection else 0

    sThreshold_area = f"{dThreshold_area:.1E}"
    sFilter = f"Lake_area > {sThreshold_area}"
    pLayer_in.SetAttributeFilter(sFilter)

    count = 0
    for pFeature_in in pLayer_in:
        pGeometry_in = pFeature_in.GetGeometryRef()
        if pGeometry_in is None:
            continue

        if iProjection == 1:
            pGeometry_in.TransformTo(srs)

        exterior_ring = pGeometry_in.GetGeometryRef(0)
        if exterior_ring is None:
            continue
        pGeometry_out = ogr.Geometry(ogr.wkbPolygon)
        pGeometry_out.AddGeometry(exterior_ring)

        aCoords_gcs = get_geometry_coordinates(pGeometry_in)
        dArea = calculate_polygon_area(aCoords_gcs[:, 0], aCoords_gcs[:, 1]) / 1.0e6

        if dArea < dThreshold_area:
            continue

        pFeature_out.SetGeometry(pGeometry_out)
        pFeature_out.SetField("area", float(dArea))
        pLayer_out.CreateFeature(pFeature_out)
        count += 1

    pDataSource_out = None
    pDataSource_in = None
    logger.info("Created %d lake polygons in %s", count, sFilename_out)
    return sFilename_out


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 2:
        create_lake_polygon(sys.argv[1], sys.argv[2])
