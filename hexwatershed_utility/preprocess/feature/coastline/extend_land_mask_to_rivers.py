import logging
import os
from typing import Sequence

from osgeo import ogr
from pyearth.gis.gdal.gdal_vector_format_support import get_vector_driver_from_filename
from pyearth.gis.geometry.calculate_distance_based_on_longitude_latitude import (
    calculate_distance_based_on_longitude_latitude,
)
from pyearth.toolbox.management.vector.merge_features import merge_features
from pyearth.toolbox.management.vector.remove_internal_polygon import remove_internal_polygon
from pyearthbuffer import gcsbuffer

try:
    from rtree.index import Index as RTreeindex

    HAS_RTREE = True
except ImportError:
    RTreeindex = None
    HAS_RTREE = False

logger = logging.getLogger(__name__)


def _bbox_overlap(bbox1, bbox2, tolerance=1e-10):
    """Check if two bounding boxes (minX, maxX, minY, maxY) overlap."""
    return not (
        bbox1[1] + tolerance < bbox2[0]
        or bbox2[1] + tolerance < bbox1[0]
        or bbox1[3] + tolerance < bbox2[2]
        or bbox2[3] + tolerance < bbox1[2]
    )


def extend_land_mask_to_rivers(
    aFilename_flowlines: Sequence[str],
    sFilename_vector_land: str,
    sFilename_vector_land_updated: str,
    iFlag_remove_internal_rings: bool = False,
) -> str:
    """
    Extend a land polygon mask to encompass coastal river flowlines that protrude into the ocean.

    Identifies flowline segments that lie outside the land polygon boundary, generates
    adaptive buffer zones around them, and unions the buffer zones with the land polygon.
    Ensures complete coverage of river mouths in land-surface and mesh models.

    Parameters
    ----------
    aFilename_flowlines : Sequence[str]
        List of river flowline vector file paths (e.g. HydroRIVERS, NHDPlus, MERIT-Hydro).
    sFilename_vector_land : str
        Path to input land polygon vector file (e.g. Natural Earth or custom coastline).
    sFilename_vector_land_updated : str
        Path to output updated land polygon vector file.
    iFlag_remove_internal_rings : bool, default=False
        If True, removes interior rings (inland holes) after merging buffer zones.

    Returns
    -------
    str
        Path to the updated land polygon file.
    """
    if not os.path.exists(sFilename_vector_land):
        raise FileNotFoundError(f"Input land vector file not found: {sFilename_vector_land}")

    pDataset_land = ogr.Open(sFilename_vector_land, 0)
    if pDataset_land is None:
        raise RuntimeError(f"Could not open land vector dataset: {sFilename_vector_land}")

    pLayer_land = pDataset_land.GetLayer(0)
    if pLayer_land is None:
        raise RuntimeError("No layer found in land vector file")

    pSpatial_reference = pLayer_land.GetSpatialRef()

    # Build spatial index for land polygons
    logger.info("Indexing land polygons...")
    land_geometries = []

    if HAS_RTREE:
        land_spatial_index = RTreeindex()
    else:
        land_spatial_index = {}

    pLayer_land.ResetReading()
    for idx, pFeature_land in enumerate(pLayer_land):
        pGeometry_land = pFeature_land.GetGeometryRef()
        if pGeometry_land is None:
            continue

        land_geometries.append(pGeometry_land.Clone())
        envelope = pGeometry_land.GetEnvelope()  # (minX, maxX, minY, maxY)

        if HAS_RTREE:
            # RTree expects (minX, minY, maxX, maxY)
            bbox = (envelope[0], envelope[2], envelope[1], envelope[3])
            land_spatial_index.insert(idx, bbox)
        else:
            land_spatial_index[idx] = envelope

    logger.info("Indexed %d land polygon features", len(land_geometries))

    aBuffer_geometries = []
    nFlowline_files = len(aFilename_flowlines)

    for iFlowline, sFilename_flowline in enumerate(aFilename_flowlines):
        if not os.path.exists(sFilename_flowline):
            logger.warning("Flowline file does not exist: %s", sFilename_flowline)
            continue

        logger.info(
            "Processing flowlines %d/%d: %s",
            iFlowline + 1,
            nFlowline_files,
            os.path.basename(sFilename_flowline),
        )
        pDataset_flowline = ogr.Open(sFilename_flowline, 0)
        if pDataset_flowline is None:
            logger.warning("Could not open flowline file: %s", sFilename_flowline)
            continue

        pLayer_flowline = pDataset_flowline.GetLayer(0)
        if pLayer_flowline is None:
            pDataset_flowline = None
            continue

        pLayer_flowline.ResetReading()
        for pFeature_flowline in pLayer_flowline:
            pGeometry_flowline = pFeature_flowline.GetGeometryRef()
            if pGeometry_flowline is None:
                continue

            flowline_envelope = pGeometry_flowline.GetEnvelope()

            candidate_indices = []
            if HAS_RTREE:
                query_bbox = (
                    flowline_envelope[0],
                    flowline_envelope[2],
                    flowline_envelope[1],
                    flowline_envelope[3],
                )
                candidate_indices = list(land_spatial_index.intersection(query_bbox))
            else:
                for idx, land_envelope in land_spatial_index.items():
                    if _bbox_overlap(flowline_envelope, land_envelope):
                        candidate_indices.append(idx)

            bIntersects_land = False
            whole_flowline = pGeometry_flowline.Clone()
            whole_flowline.FlattenTo2D()

            for candidate_idx in candidate_indices:
                if candidate_idx >= len(land_geometries):
                    continue
                if whole_flowline is None or whole_flowline.IsEmpty():
                    break

                pGeometry_land = land_geometries[candidate_idx]
                if pGeometry_land is None:
                    continue

                if whole_flowline.Intersects(pGeometry_land):
                    bIntersects_land = True
                    outside_part = whole_flowline.Difference(pGeometry_land)
                    if outside_part is not None and not outside_part.IsEmpty():
                        outside_part.FlattenTo2D()
                        whole_flowline = outside_part
                    else:
                        whole_flowline = None
                        break

            lines_to_buffer = []
            if bIntersects_land and whole_flowline is not None and not whole_flowline.IsEmpty():
                outside_flowline = whole_flowline
                sGeometry_type = outside_flowline.GetGeometryName().upper()
                if sGeometry_type == "LINESTRING":
                    lines_to_buffer.append(outside_flowline)
                elif sGeometry_type == "MULTILINESTRING":
                    for k in range(outside_flowline.GetGeometryCount()):
                        sub_line = outside_flowline.GetGeometryRef(k)
                        if sub_line is not None and not sub_line.IsEmpty():
                            lines_to_buffer.append(sub_line)
            elif not bIntersects_land:
                outside_flowline = pGeometry_flowline.Clone()
                outside_flowline.FlattenTo2D()
                sGeometry_type = outside_flowline.GetGeometryName().upper()
                if sGeometry_type == "LINESTRING":
                    lines_to_buffer.append(outside_flowline)
                elif sGeometry_type == "MULTILINESTRING":
                    for k in range(outside_flowline.GetGeometryCount()):
                        sub_line = outside_flowline.GetGeometryRef(k)
                        if sub_line is not None and not sub_line.IsEmpty():
                            lines_to_buffer.append(sub_line)

            for line in lines_to_buffer:
                pt_count = line.GetPointCount()
                if pt_count < 2:
                    continue
                p_start = line.GetPoint(0)
                p_end = line.GetPoint(pt_count - 1)
                distance = calculate_distance_based_on_longitude_latitude(
                    p_start[0], p_start[1], p_end[0], p_end[1]
                )
                if distance <= 0:
                    continue
                buf_distance = max(distance / 2.0, 1.0)
                pBuf_geom = gcsbuffer(line, distance=buf_distance, output_format="ogr")
                if pBuf_geom is not None and not pBuf_geom.IsEmpty():
                    aBuffer_geometries.append(pBuf_geom)

        pDataset_flowline = None

    # Write temporary file combining original features + buffer geometries
    sDir_out = os.path.dirname(sFilename_vector_land_updated)
    if sDir_out and not os.path.exists(sDir_out):
        os.makedirs(sDir_out, exist_ok=True)

    sExt = os.path.splitext(sFilename_vector_land_updated)[1] or ".geojson"
    sBase = os.path.splitext(os.path.basename(sFilename_vector_land_updated))[0]
    sTmp_path = os.path.join(sDir_out or ".", f"{sBase}_tmp{sExt}")

    pDriver = get_vector_driver_from_filename(sTmp_path)
    if pDriver is None:
        pDriver = ogr.GetDriverByName("GeoJSON")

    if os.path.exists(sTmp_path):
        pDriver.DeleteDataSource(sTmp_path)

    pDataset_output = pDriver.CreateDataSource(sTmp_path)
    pLayer_output = pDataset_output.CreateLayer(
        "land_updated", pSpatial_reference, geom_type=ogr.wkbPolygon
    )

    # Copy fields
    pLayerDefn_input = pLayer_land.GetLayerDefn()
    for i in range(pLayerDefn_input.GetFieldCount()):
        pLayer_output.CreateField(pLayerDefn_input.GetFieldDefn(i))

    pLayer_land.ResetReading()
    for pFeature in pLayer_land:
        pFeature_out = ogr.Feature(pLayer_output.GetLayerDefn())
        pFeature_out.SetFrom(pFeature)
        pLayer_output.CreateFeature(pFeature_out)
        pFeature_out = None

    # Append coastal buffer polygons
    for pBuffer_geom in aBuffer_geometries:
        pFeature_out = ogr.Feature(pLayer_output.GetLayerDefn())
        pFeature_out.SetGeometry(pBuffer_geom)
        pLayer_output.CreateFeature(pFeature_out)
        pFeature_out = None

    pDataset_output = None
    pDataset_land = None

    # Merge / dissolve overlapping features into clean output
    if iFlag_remove_internal_rings:
        sTmp_merged = os.path.join(sDir_out or ".", f"{sBase}_tmp_merged{sExt}")
        merge_features(sTmp_path, sTmp_merged)
        remove_internal_polygon(sTmp_merged, sFilename_vector_land_updated)
        if os.path.exists(sTmp_merged):
            try:
                pDriver.DeleteDataSource(sTmp_merged)
            except Exception:
                pass
    else:
        merge_features(sTmp_path, sFilename_vector_land_updated)

    if os.path.exists(sTmp_path):
        try:
            pDriver.DeleteDataSource(sTmp_path)
        except Exception:
            pass

    logger.info(
        "Added %d coastal buffer zones. Updated land mask: %s",
        len(aBuffer_geometries),
        sFilename_vector_land_updated,
    )
    return sFilename_vector_land_updated


def geometries_bbox_overlap(bbox1: tuple, bbox2: tuple, tolerance: float = 1e-10) -> bool:
    """Check if two bounding boxes (minX, maxX, minY, maxY) overlap with optional tolerance."""
    return _bbox_overlap(bbox1, bbox2, tolerance=tolerance)


__all__ = [
    "extend_land_mask_to_rivers",
    "geometries_bbox_overlap",
]
