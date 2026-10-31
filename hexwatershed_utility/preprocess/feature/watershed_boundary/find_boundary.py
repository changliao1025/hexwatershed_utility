import glob
import logging
import os
from typing import List, Optional, Sequence, Tuple, Union

from osgeo import ogr
from pyearth.gis.gdal.read.vector.gdal_get_vector_extent import gdal_get_vector_extent
from pyearth.gis.location.get_hydrosheds_continent_from_extent import (
    get_hydrosheds_continent_from_extent,
)

logger = logging.getLogger(__name__)


def _find_containing_feature_in_file(
    pGeometry_target: ogr.Geometry, sFilename_boundary: str
) -> Optional[str]:
    """
    Search a boundary vector file for a polygon that contains pGeometry_target.

    Args:
        pGeometry_target (ogr.Geometry): Target geometry to search for.
        sFilename_boundary (str): Path to boundary vector file (Shapefile, GeoJSON, etc.).

    Returns:
        Optional[str]:
            - WKT string if a containing polygon is found.
            - "INTERSECTS" if polygons intersect but none completely contain the target.
            - None if neither containment nor intersection is found.
    """
    try:
        pDataset_boundary = ogr.Open(sFilename_boundary, 0)
        if pDataset_boundary is None:
            logger.error("Could not open boundary file: %s", sFilename_boundary)
            return None

        pLayer_boundary = pDataset_boundary.GetLayer(0)
        if pLayer_boundary is None:
            logger.error("Could not get layer from boundary file: %s", sFilename_boundary)
            pDataset_boundary = None
            return None

        # Apply spatial bounding box filter if possible for performance
        try:
            envelope = pGeometry_target.GetEnvelope()  # (minX, maxX, minY, maxY)
            pLayer_boundary.SetSpatialFilterRect(envelope[0], envelope[2], envelope[1], envelope[3])
        except Exception as e:
            logger.debug("Could not set spatial filter rect: %s", e)
            pLayer_boundary.SetSpatialFilter(None)

        bIntersects_any = False
        pLayer_boundary.ResetReading()
        for pFeature_boundary in pLayer_boundary:
            pGeometry_boundary = pFeature_boundary.GetGeometryRef()
            if pGeometry_boundary is None:
                continue

            try:
                if pGeometry_boundary.Contains(pGeometry_target):
                    wkt = pGeometry_boundary.ExportToWkt()
                    pDataset_boundary = None
                    return wkt
                elif pGeometry_boundary.Intersects(pGeometry_target):
                    bIntersects_any = True
            except Exception as e:
                logger.debug("Error checking geometry containment/intersection: %s", e)
                continue

        pDataset_boundary = None
        if bIntersects_any:
            return "INTERSECTS"
        return None

    except Exception as e:
        logger.error("Error processing boundary file %s: %s", sFilename_boundary, e)
        return None


def _find_boundary_for_geometry(
    pGeometry_target: ogr.Geometry,
    aBoundary_files: Sequence[str],
) -> Optional[str]:
    """
    Search an ordered sequence of boundary files (from fine to coarse scale)
    to find the smallest containing polygon for pGeometry_target.
    """
    for sFilename in aBoundary_files:
        if not os.path.exists(sFilename):
            logger.warning("Boundary file does not exist: %s", sFilename)
            continue

        sResult = _find_containing_feature_in_file(pGeometry_target, sFilename)
        if sResult == "INTERSECTS":
            # Target geometry spans across boundaries at this scale; ascend to coarser scale
            continue
        elif sResult is not None:
            logger.info("Found containing boundary in: %s", sFilename)
            return sResult

    logger.warning("No containing boundary found across provided files")
    return None


def _discover_boundary_files(
    sFolder_boundary: str,
    sDataset_type: str = "hydrobasins",
    aExtent: Optional[Tuple[float, float, float, float]] = None,
    aLevels: Optional[Sequence[Union[int, str]]] = None,
    sPattern: Optional[str] = None,
) -> List[str]:
    """
    Discover and order boundary files from fine scale to coarse scale.

    Args:
        sFolder_boundary (str): Base directory containing boundary files.
        sDataset_type (str): Type of dataset ("hydrobasins", "generic", or "custom").
        aExtent (Optional[Tuple]): Vector extent (minx, maxx, miny, maxy) for regional resolution.
        aLevels (Optional[Sequence]): Scale/hierarchy levels ordered from fine to coarse.
        sPattern (Optional[str]): Custom glob pattern.

    Returns:
        List[str]: Ordered list of boundary file paths.
    """
    aFiles = []

    if sDataset_type.lower() in ("hydrobasins", "hydrosheds"):
        # Resolve regional/continent directory
        sFolder_search = sFolder_boundary
        if aExtent is not None:
            sRegion = get_hydrosheds_continent_from_extent(aExtent)
            aMatched_folders = glob.glob(os.path.join(sFolder_boundary, f"*_{sRegion}_*"))
            for sPath in aMatched_folders:
                if os.path.isdir(sPath):
                    sFolder_search = sPath
                    break

        levels: Sequence[Union[int, str]]
        if aLevels is None:
            levels = list(range(10, 0, -1))
        else:
            levels = aLevels

        for level in levels:
            sLevel = f"lev{int(level):02d}" if isinstance(level, int) else f"lev{level}"
            pattern = f"hybas_lake_*{sLevel}*_*.shp"
            matches = glob.glob(os.path.join(sFolder_search, pattern))
            if matches:
                aFiles.append(matches[0])
            else:
                logger.debug("No shapefile found for level %s in %s", sLevel, sFolder_search)

    else:
        # Generic or custom dataset
        if sPattern is not None:
            aFiles = sorted(glob.glob(os.path.join(sFolder_boundary, sPattern)))
        elif aLevels is not None:
            for custom_level in aLevels:
                pattern = f"*{custom_level}*"
                matches = sorted(glob.glob(os.path.join(sFolder_boundary, pattern)))
                aFiles.extend(matches)
        else:
            for ext in ("*.shp", "*.geojson", "*.gpkg", "*.parquet"):
                aFiles.extend(sorted(glob.glob(os.path.join(sFolder_boundary, ext))))

    return aFiles


def find_containing_boundary(
    sFilename_river_network_in: str,
    sFolder_boundary_in: Optional[str] = None,
    aBoundary_files: Optional[Sequence[str]] = None,
    sDataset_type: str = "hydrobasins",
    iFlag_nested_in: bool = False,
    aLevels: Optional[Sequence[Union[int, str]]] = None,
    sPattern: Optional[str] = None,
) -> Optional[Union[str, List[str]]]:
    """
    Find the minimal/smallest watershed boundary that completely contains a river network
    or its constituent subbasins from hierarchical boundary datasets.

    Supports global datasets (e.g. HydroBASINS) as well as regional or custom datasets
    (e.g. MERIT-Hydro, NHDPlus, or user-supplied boundary files).

    Args:
        sFilename_river_network_in (str): Path to input river network vector file (e.g. GeoJSON).
        sFolder_boundary_in (Optional[str]): Path to folder containing hierarchical boundary files.
        aBoundary_files (Optional[Sequence[str]]): Direct list of boundary file paths ordered
            from fine scale to coarse scale. If provided, overrides folder discovery.
        sDataset_type (str): Dataset provider/convention ("hydrobasins", "generic", "custom").
        iFlag_nested_in (bool): If True, process each subbasin feature separately and return
            a list of containing boundary WKTs. If False, process the union of all features
            and return a single containing boundary WKT.
        aLevels (Optional[Sequence]): Hierarchy levels to search (from fine to coarse).
        sPattern (Optional[str]): Custom glob pattern for boundary file discovery.

    Returns:
        Optional[Union[str, List[str]]]:
            - If iFlag_nested_in is False: WKT string of containing boundary, or None.
            - If iFlag_nested_in is True: List of WKT strings for each subbasin, or None.
    """
    logger.info("Reading river network from: %s", sFilename_river_network_in)

    # Read vector extent and union geometry
    aNetwork_bounds, pGeometry_network = gdal_get_vector_extent(
        sFilename_river_network_in, iFlag_return_union_geometry=1
    )
    if aNetwork_bounds is None:
        logger.error("Failed to read extent from: %s", sFilename_river_network_in)
        return None

    # Resolve boundary files
    if aBoundary_files is not None:
        files_to_search = list(aBoundary_files)
    elif sFolder_boundary_in is not None:
        files_to_search = _discover_boundary_files(
            sFolder_boundary=sFolder_boundary_in,
            sDataset_type=sDataset_type,
            aExtent=aNetwork_bounds,
            aLevels=aLevels,
            sPattern=sPattern,
        )
    else:
        raise ValueError("Either `aBoundary_files` or `sFolder_boundary_in` must be provided.")

    if not files_to_search:
        logger.error("No boundary files found to search.")
        return None

    logger.info("Searching through %d boundary file(s)", len(files_to_search))

    if iFlag_nested_in:
        # Process each feature/subbasin individually
        aBoundaries = []
        try:
            pDataset_network = ogr.Open(sFilename_river_network_in, 0)
            if pDataset_network is None:
                logger.error("Could not open river network file: %s", sFilename_river_network_in)
                return None

            pLayer_network = pDataset_network.GetLayer(0)
            if pLayer_network is None:
                logger.error(
                    "Could not get layer from river network: %s", sFilename_river_network_in
                )
                return None

            pLayer_network.ResetReading()
            for pFeature_network in pLayer_network:
                pGeometry_feature = pFeature_network.GetGeometryRef()
                if pGeometry_feature is None:
                    continue

                sWkt = _find_boundary_for_geometry(pGeometry_feature, files_to_search)
                if sWkt:
                    aBoundaries.append(sWkt)
                    logger.info("Boundary found for feature FID %s", pFeature_network.GetFID())

            pDataset_network = None

        except Exception as e:
            logger.error("Error processing nested features: %s", e)
            return None

        if not aBoundaries:
            logger.warning("No containing boundary found for any subbasin feature")
            return None

        return aBoundaries

    else:
        # Process union geometry
        return _find_boundary_for_geometry(pGeometry_network, files_to_search)


__all__ = [
    "find_containing_boundary",
]
