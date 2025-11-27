from osgeo import gdal, ogr, osr
from pyearth.gis.gdal.gdal_vector_format_support import get_vector_driver_from_filename

def reorganize_boundary(aWkt_in, sFilename_out=None):
    """
    Reorganize the watershed boundary WKT by handling overlapping polygons and extracting
    minimal polygon components through intersection and difference operations.

    Args:
        aWkt_in (list): List of WKT strings representing the watershed boundary polygons.

    Returns:
        sWkt_out (str): WKT string of the reorganized watershed boundary multipolygon
                       containing all minimal polygon components.
    """

    # Create a spatial reference for WGS84
    oSRS = osr.SpatialReference()
    oSRS.ImportFromEPSG(4326)

    # Parse input geometries and validate
    aPolygons = []
    for sWkt in aWkt_in:
        oGeom = ogr.CreateGeometryFromWkt(sWkt)
        if oGeom is None:
            continue

        # Ensure the geometry is valid and a polygon
        if oGeom.GetGeometryType() == ogr.wkbPolygon:
            # Ensure proper orientation (exterior ring clockwise, holes counter-clockwise)

            aPolygons.append(oGeom.Clone())
        elif oGeom.GetGeometryType() == ogr.wkbMultiPolygon:
            for i in range(oGeom.GetGeometryCount()):
                oSubGeom = oGeom.GetGeometryRef(i)
                if oSubGeom.GetGeometryType() == ogr.wkbPolygon:
                    aPolygons.append(oSubGeom.Clone())

    if not aPolygons:
        # Return empty multipolygon if no valid polygons
        oFinalGeom = ogr.Geometry(ogr.wkbMultiPolygon)
        return oFinalGeom.ExportToWkt()

    # Extract minimal polygons by decomposing overlaps
    aMinimalPolygons = _extract_minimal_polygons(aPolygons)

    # Create final multipolygon from minimal components
    oFinalGeom = ogr.Geometry(ogr.wkbMultiPolygon)
    for oPolygon in aMinimalPolygons:
        if oPolygon.IsValid() and not oPolygon.IsEmpty():
            oFinalGeom.AddGeometry(oPolygon)

    # Return the WKT of the reorganized boundary
    sWkt_out = oFinalGeom.ExportToWkt()

    if sFilename_out is not None:
        # Save to file if output filename is provided
        pDriver = get_vector_driver_from_filename(sFilename_out)
        if pDriver is None:
            raise ValueError(f"Could not get driver for file: {sFilename_out}")

        pDataset_out = pDriver.CreateDataSource(sFilename_out)
        if pDataset_out is None:
            raise ValueError(f"Could not create output dataset: {sFilename_out}")

        pLayer_out = pDataset_out.CreateLayer("watershed_boundary", geom_type=ogr.wkbMultiPolygon, srs=oSRS)
        if pLayer_out is None:
            raise ValueError(f"Could not create layer in output dataset: {sFilename_out}")

        pFeature_out = ogr.Feature(pLayer_out.GetLayerDefn())
        pGeometry_out = ogr.CreateGeometryFromWkt(sWkt_out)
        if pGeometry_out is None:
            raise ValueError(f"Could not create geometry from WKT: {sWkt_out}")

        pFeature_out.SetGeometry(pGeometry_out)
        if pLayer_out.CreateFeature(pFeature_out) != ogr.OGRERR_NONE:
            raise ValueError(f"Could not create feature in output layer: {sFilename_out}")

        # Clean up
        pFeature_out = None
        pDataset_out = None
    return sWkt_out


def _extract_minimal_polygons(aPolygons):
    """
    Remove duplicate polygons and keep only unique ones.

    This function filters out duplicate polygons based on geometric equality,
    keeping only unique polygon geometries.

    Args:
        aPolygons (list): List of OGR Geometry objects (polygons).

    Returns:
        aUniquePolygons (list): List of unique polygon components.
    """
    if len(aPolygons) <= 1:
        return aPolygons

    aUniquePolygons = []

    for i, oCurrentPoly in enumerate(aPolygons):
        if oCurrentPoly.IsEmpty():
            continue

        bIsDuplicate = False

        # Check if this polygon is a duplicate of any already processed polygon
        for oUniquePoly in aUniquePolygons:
            try:
                # Check if polygons are geometrically equal
                if oCurrentPoly.Equals(oUniquePoly):
                    bIsDuplicate = True
                    break

                # Also check if they have the same area and centroid (for numerical precision issues)
                if (abs(oCurrentPoly.GetArea() - oUniquePoly.GetArea()) < 1e-10 and
                    oCurrentPoly.Centroid().Equals(oUniquePoly.Centroid())):
                    bIsDuplicate = True
                    break

            except Exception as e:
                print(f"Error comparing polygons: {e}")
                continue

        # Add polygon if it's not a duplicate
        if not bIsDuplicate:
            aUniquePolygons.append(oCurrentPoly.Clone())

    return aUniquePolygons

