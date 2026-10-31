import logging
import os
from typing import Optional, Sequence

from osgeo import ogr, osr
from pyearth.gis.gdal.gdal_vector_format_support import get_vector_driver_from_filename
from pyearth.toolbox.management.vector.merge_features import merge_features
from pyearth.toolbox.management.vector.merge_vector_files import merge_vector_files

logger = logging.getLogger(__name__)


def _extract_unique_polygons(aPolygons):
    """
    Remove duplicate polygons and keep only unique polygon geometries.

    Uses envelope and area pre-filters followed by geometric equality checks.

    Args:
        aPolygons (list): List of OGR Geometry objects (polygons).

    Returns:
        list: List of unique polygon components.
    """
    if len(aPolygons) <= 1:
        return aPolygons

    aUniquePolygons = []

    for oCurrentPoly in aPolygons:
        if oCurrentPoly.IsEmpty():
            continue

        bIsDuplicate = False
        env_curr = oCurrentPoly.GetEnvelope()
        area_curr = oCurrentPoly.GetArea()

        for oUniquePoly in aUniquePolygons:
            env_u = oUniquePoly.GetEnvelope()
            # Fast envelope rejection
            if (
                abs(env_curr[0] - env_u[0]) > 1e-7
                or abs(env_curr[1] - env_u[1]) > 1e-7
                or abs(env_curr[2] - env_u[2]) > 1e-7
                or abs(env_curr[3] - env_u[3]) > 1e-7
            ):
                continue

            # Fast area rejection
            area_u = oUniquePoly.GetArea()
            if abs(area_curr - area_u) > 1e-7:
                continue

            try:
                # Geometric equality check
                if oCurrentPoly.Equals(oUniquePoly):
                    bIsDuplicate = True
                    break

                # Centroid check for precision issues
                c1 = oCurrentPoly.Centroid()
                c2 = oUniquePoly.Centroid()
                if abs(c1.GetX() - c2.GetX()) < 1e-7 and abs(c1.GetY() - c2.GetY()) < 1e-7:
                    bIsDuplicate = True
                    break

            except Exception as e:
                logger.debug("Error comparing polygons: %s", e)
                continue

        if not bIsDuplicate:
            aUniquePolygons.append(oCurrentPoly.Clone())

    return aUniquePolygons


def merge_boundary_polygons(
    aWkt_in: Sequence[str],
    sFilename_out: Optional[str] = None,
    iFlag_union: bool = True,
) -> str:
    """
    Deduplicate and merge watershed boundary polygon WKTs into a single MultiPolygon,
    optionally unioning contiguous/overlapping polygons into a single dissolved feature.

    Parameters
    ----------
    aWkt_in : Sequence[str]
        List of WKT strings representing boundary polygons.
    sFilename_out : Optional[str], default=None
        Optional output vector file path to save the merged result.
    iFlag_union : bool, default=True
        If True, union and dissolve overlapping or contiguous polygons into a single
        unified boundary. If False, aggregate deduplicated polygons without dissolving.

    Returns
    -------
    str
        WKT string of the merged MultiPolygon containing unique polygon components.
    """
    oSRS = osr.SpatialReference()
    oSRS.ImportFromEPSG(4326)

    # Parse input geometries and validate
    aPolygons = []
    for sWkt in aWkt_in:
        oGeom = ogr.CreateGeometryFromWkt(sWkt)
        if oGeom is None:
            continue

        if oGeom.GetGeometryType() == ogr.wkbPolygon:
            aPolygons.append(oGeom.Clone())
        elif oGeom.GetGeometryType() == ogr.wkbMultiPolygon:
            for i in range(oGeom.GetGeometryCount()):
                oSubGeom = oGeom.GetGeometryRef(i)
                if oSubGeom.GetGeometryType() == ogr.wkbPolygon:
                    aPolygons.append(oSubGeom.Clone())

    if not aPolygons:
        oFinalGeom = ogr.Geometry(ogr.wkbMultiPolygon)
        return oFinalGeom.ExportToWkt()

    # Deduplicate unique polygons
    aUniquePolygons = _extract_unique_polygons(aPolygons)

    if iFlag_union and len(aUniquePolygons) > 1:
        oUnionGeom = aUniquePolygons[0].Clone()
        for oPoly in aUniquePolygons[1:]:
            try:
                oMerged = oUnionGeom.Union(oPoly)
                if oMerged is not None and not oMerged.IsEmpty():
                    oUnionGeom = oMerged
            except Exception as e:
                logger.debug("Union operation failed between polygons: %s", e)

        if oUnionGeom.GetGeometryType() == ogr.wkbPolygon:
            oFinalGeom = ogr.Geometry(ogr.wkbMultiPolygon)
            oFinalGeom.AddGeometry(oUnionGeom)
        elif oUnionGeom.GetGeometryType() == ogr.wkbMultiPolygon:
            oFinalGeom = oUnionGeom
        else:
            oFinalGeom = ogr.Geometry(ogr.wkbMultiPolygon)
            for oPolygon in aUniquePolygons:
                if oPolygon.IsValid() and not oPolygon.IsEmpty():
                    oFinalGeom.AddGeometry(oPolygon)
    else:
        oFinalGeom = ogr.Geometry(ogr.wkbMultiPolygon)
        for oPolygon in aUniquePolygons:
            if oPolygon.IsValid() and not oPolygon.IsEmpty():
                oFinalGeom.AddGeometry(oPolygon)

    sWkt_out = oFinalGeom.ExportToWkt()

    if sFilename_out is not None:
        pDriver = get_vector_driver_from_filename(sFilename_out)
        if pDriver is None:
            raise ValueError(f"Could not get driver for file: {sFilename_out}")

        pDataset_out = pDriver.CreateDataSource(sFilename_out)
        if pDataset_out is None:
            raise ValueError(f"Could not create output dataset: {sFilename_out}")

        pLayer_out = pDataset_out.CreateLayer(
            "watershed_boundary", geom_type=ogr.wkbMultiPolygon, srs=oSRS
        )
        if pLayer_out is None:
            raise ValueError(f"Could not create layer in output dataset: {sFilename_out}")

        pFeature_out = ogr.Feature(pLayer_out.GetLayerDefn())
        pGeometry_out = ogr.CreateGeometryFromWkt(sWkt_out)
        if pGeometry_out is None:
            raise ValueError(f"Could not create geometry from WKT: {sWkt_out}")

        pFeature_out.SetGeometry(pGeometry_out)
        if pLayer_out.CreateFeature(pFeature_out) != ogr.OGRERR_NONE:
            raise ValueError(f"Could not create feature in output layer: {sFilename_out}")

        pFeature_out = None
        pDataset_out = None

    return sWkt_out


def merge_boundary_files(
    aBoundary_files: Sequence[str],
    sFilename_out: str,
    iFlag_union: bool = True,
) -> str:
    """
    Merge multiple boundary polygon vector files into a single output file,
    unioning/dissolving features into a single continuous watershed boundary.

    Uses `pyearth`'s `merge_vector_files` and `merge_features` APIs to union
    and dissolve internal boundaries between subbasins.

    Parameters
    ----------
    aBoundary_files : Sequence[str]
        List of input polygon vector file paths.
    sFilename_out : str
        Output file path for merged boundary (e.g. GeoJSON, Shapefile, Parquet).
    iFlag_union : bool, default=True
        If True, union and dissolve all subbasin features into a single continuous
        watershed boundary feature using pyearth.merge_features.
        If False, concatenate all features into one file without dissolving.

    Returns
    -------
    str
        Path to output merged vector file.
    """
    if not aBoundary_files:
        raise ValueError("Input aBoundary_files list cannot be empty")

    sDir_out = os.path.dirname(sFilename_out)
    if sDir_out and not os.path.exists(sDir_out):
        os.makedirs(sDir_out, exist_ok=True)

    if not iFlag_union:
        merge_vector_files(aBoundary_files, sFilename_out)
        return sFilename_out

    if len(aBoundary_files) == 1:
        merge_features(aBoundary_files[0], sFilename_out)
        return sFilename_out

    # Multiple files: merge into intermediate file, then union features into one
    sExt = os.path.splitext(sFilename_out)[1] or ".geojson"
    sBase = os.path.splitext(os.path.basename(sFilename_out))[0]
    sFilename_tmp = os.path.join(sDir_out or ".", f"{sBase}_raw_merge{sExt}")

    try:
        merge_vector_files(aBoundary_files, sFilename_tmp)
        merge_features(sFilename_tmp, sFilename_out)
    finally:
        if os.path.exists(sFilename_tmp):
            pDriver = get_vector_driver_from_filename(sFilename_tmp)
            if pDriver is not None:
                try:
                    pDriver.DeleteDataSource(sFilename_tmp)
                except Exception:
                    if os.path.exists(sFilename_tmp):
                        os.remove(sFilename_tmp)
            elif os.path.exists(sFilename_tmp):
                os.remove(sFilename_tmp)

    logger.info(
        "Merged and unioned %d boundary file(s) into: %s", len(aBoundary_files), sFilename_out
    )
    return sFilename_out


__all__ = [
    "merge_boundary_polygons",
    "merge_boundary_files",
]
