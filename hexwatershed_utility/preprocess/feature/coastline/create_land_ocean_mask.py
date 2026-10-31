import glob
import logging
import os
from typing import List, Optional, Sequence, Tuple

from pyearth.toolbox.analysis.image.raster_process import (
    create_raster_buffer_zone,
    fix_raster_antimeridian_issue,
)
from pyearth.toolbox.conversion.convert_vector_to_global_raster import (
    convert_vector_to_global_raster,
)
from pyearth.toolbox.data.ocean.define_land_ocean_mask import (
    create_land_ocean_vector_mask_naturalearth,
)
from pyearth.toolbox.management.vector.merge_features import merge_features
from pyearth.toolbox.management.vector.merge_vector_files import merge_vector_files
from pyearth.toolbox.management.vector.remove_internal_polygon import remove_internal_polygon
from pyearth.toolbox.management.vector.remove_small_polygon import remove_small_polygon
from pyearthbuffer import gcsbuffer_file

logger = logging.getLogger(__name__)


def create_land_ocean_mask_from_basins(
    sWorkspace_output: str,
    aBoundary_files: Optional[Sequence[str]] = None,
    sFolder_boundary_in: Optional[str] = None,
    sLevel: str = "lev03",
    dResolution_x: float = 30 / 3600.0,
    dResolution_y: float = 30 / 3600.0,
    dThreshold_area_island: float = 1000.0,
    dResolution_coastline_buffer: float = 10.0,
) -> Tuple[str, str]:
    """
    Create a global land-ocean mask from regional / continental watershed basin boundaries.

    Merges regional/continental basin boundaries, removes small islands,
    generates a coastal buffer zone, and rasterizes to a global GeoTIFF.

    Parameters
    ----------
    sWorkspace_output : str
        Directory to save intermediate and final mask outputs.
    aBoundary_files : Optional[Sequence[str]], default=None
        Direct list of boundary vector files to merge. If provided, overrides folder discovery.
    sFolder_boundary_in : Optional[str], default=None
        Directory containing hierarchical boundary datasets (e.g. HydroBASINS).
    sLevel : str, default='lev03'
        Hierarchy level tag for folder discovery (e.g. 'lev03').
    dResolution_x : float, default=30/3600.0
        X-resolution in degrees.
    dResolution_y : float, default=30/3600.0
        Y-resolution in degrees.
    dThreshold_area_island : float, default=1000.0
        Minimum area threshold in km²; smaller disconnected polygons/islands are filtered out.
    dResolution_coastline_buffer : float, default=10.0
        Coastline buffer distance in km.

    Returns
    -------
    Tuple[str, str]
        (Path to output raster mask GeoTIFF, Path to output merged vector Parquet).
    """
    os.makedirs(sWorkspace_output, exist_ok=True)

    # Determine input boundary files
    if aBoundary_files is not None:
        boundary_files = list(aBoundary_files)
    elif sFolder_boundary_in is not None:
        boundary_files = []
        # Support HydroBASINS or general regional folder discovery
        regions = ["af", "ar", "as", "au", "eu", "gr", "na", "sa", "si"]
        for sRegion in regions:
            sFolder = f"hybas_lake_{sRegion}_lev01-12_v1c"
            sFolder_in = os.path.join(sFolder_boundary_in, sFolder)
            if not os.path.exists(sFolder_in):
                sFolder_in = sFolder_boundary_in

            pattern = f"hybas_lake_*{sLevel}*_*.shp"
            matches = glob.glob(os.path.join(sFolder_in, pattern))
            if matches:
                boundary_files.append(matches[0])

        if not boundary_files:
            # Fallback: search for any matching level files directly in folder
            boundary_files = glob.glob(os.path.join(sFolder_boundary_in, f"*{sLevel}*.shp"))
    else:
        raise ValueError("Either `aBoundary_files` or `sFolder_boundary_in` must be provided.")

    if not boundary_files:
        raise FileNotFoundError("No boundary files found for land-ocean mask generation.")

    logger.info("Merging %d boundary file(s)...", len(boundary_files))
    sFilename_merge = os.path.join(sWorkspace_output, f"basins_merged_{sLevel}.parquet")
    merge_vector_files(boundary_files, sFilename_merge)

    # Remove small island polygons
    sFilename_wo_island = os.path.join(sWorkspace_output, "land_ocean_mask_wo_island.parquet")
    remove_small_polygon(sFilename_merge, sFilename_wo_island, dThreshold_area_island)

    # Dissolve / merge features into unified land polygon
    sFilename_vector_merged_raw = os.path.join(
        sWorkspace_output, "land_ocean_mask_wo_island_merged_raw.parquet"
    )
    merge_features(sFilename_wo_island, sFilename_vector_merged_raw)

    # Rasterize to global raster GeoTIFF
    sFilename_tif_wo_island = os.path.join(sWorkspace_output, "land_ocean_mask_wo_island.tif")
    convert_vector_to_global_raster(
        sFilename_vector_merged_raw,
        sFilename_tif_wo_island,
        dResolution_x,
        dResolution_y,
        iFlag_boundary_only_in=0,
        dFill_value_in=2,
    )

    # Generate coastal buffer zone
    sFilename_parquet_buffer = os.path.join(sWorkspace_output, "land_ocean_mask_buffer.parquet")
    gcsbuffer_file(
        input_file=sFilename_vector_merged_raw,
        output_file=sFilename_parquet_buffer,
        distance=dResolution_coastline_buffer * 1000.0,
    )

    sFilename_vector_parquet_merged = os.path.join(
        sWorkspace_output, "land_ocean_mask_buffer_merged.parquet"
    )
    merge_features(sFilename_parquet_buffer, sFilename_vector_parquet_merged)

    logger.info("Completed land-ocean mask from basins: %s", sFilename_tif_wo_island)
    return sFilename_tif_wo_island, sFilename_vector_parquet_merged


def create_land_ocean_mask_from_shoreline(
    sWorkspace_output: str,
    sFilename_land_vector_in: Optional[str] = None,
    dResolution_x: float = 30 / 3600.0,
    dResolution_y: float = 30 / 3600.0,
    dThreshold_area_island: float = 1000.0,
    iRaster_buffer_pixel: int = 1,
    iFlag_remove_internal_rings: bool = False,
) -> Tuple[str, str]:
    """
    Create land-ocean mask from global shoreline vector data (e.g. Natural Earth).

    Filters small islands, optionally removes internal lakes/holes, rasterizes, and fixes
    antimeridian dateline seam issues.

    Parameters
    ----------
    sWorkspace_output : str
        Directory to save outputs.
    sFilename_land_vector_in : Optional[str], default=None
        Path to land shoreline vector file. If None, auto-downloads Natural Earth 10m land.
    dResolution_x : float, default=30/3600.0
        X-resolution in degrees.
    dResolution_y : float, default=30/3600.0
        Y-resolution in degrees.
    dThreshold_area_island : float, default=1000.0
        Minimum area threshold in km² for island filtering.
    iRaster_buffer_pixel : int, default=1
        Number of pixels to buffer the raster coastline.
    iFlag_remove_internal_rings : bool, default=False
        If True, removes interior rings (inland lakes) from land polygons so entire continents
        are treated as continuous land. If False (default), retains interior rings.

    Returns
    -------
    Tuple[str, str]
        (Path to output raster GeoTIFF with antimeridian fixed, Path to cleaned vector GeoJSON).
    """
    os.makedirs(sWorkspace_output, exist_ok=True)

    if sFilename_land_vector_in is None:
        sFilename_raw = os.path.join(sWorkspace_output, "land_ocean_mask_naturalearth.geojson")
        create_land_ocean_vector_mask_naturalearth(sFilename_raw)
    else:
        sFilename_raw = sFilename_land_vector_in

    # Filter small islands
    sFilename_wo_island = os.path.join(sWorkspace_output, "land_ocean_mask_wo_island.geojson")
    remove_small_polygon(sFilename_raw, sFilename_wo_island, dThreshold_area_island)

    # Remove internal rings (inland lakes) if requested
    if iFlag_remove_internal_rings:
        sFilename_vector_processed = os.path.join(
            sWorkspace_output, "land_ocean_mask_wo_internal_ring.geojson"
        )
        remove_internal_polygon(sFilename_wo_island, sFilename_vector_processed)
    else:
        sFilename_vector_processed = sFilename_wo_island

    # Convert to global raster
    sFilename_tif_wo_island = os.path.join(sWorkspace_output, "land_ocean_mask_wo_island.tif")
    convert_vector_to_global_raster(
        sFilename_vector_processed,
        sFilename_tif_wo_island,
        dResolution_x,
        dResolution_y,
        iFlag_boundary_only_in=0,
        dFill_value_in=2,
    )

    # Buffer raster and fix antimeridian wrapping
    if iRaster_buffer_pixel > 0:
        sFilename_buffered = os.path.join(
            sWorkspace_output, "land_ocean_mask_wo_island_buffered.tif"
        )
        create_raster_buffer_zone(
            sFilename_tif_wo_island,
            sFilename_buffered,
            1,
            int(iRaster_buffer_pixel),
        )
        sFilename_fixed = os.path.join(
            sWorkspace_output, "land_ocean_mask_wo_island_buffered_fixed.tif"
        )
        fix_raster_antimeridian_issue(
            sFilename_buffered,
            sFilename_fixed,
            1,
            2,
            iRaster_buffer_pixel,
        )
    else:
        sFilename_fixed = os.path.join(
            sWorkspace_output, "land_ocean_mask_wo_island_buffered_fixed.tif"
        )
        fix_raster_antimeridian_issue(
            sFilename_tif_wo_island,
            sFilename_fixed,
            1,
            2,
            1,
        )

    logger.info("Completed land-ocean mask from shoreline: %s", sFilename_fixed)
    return sFilename_fixed, sFilename_vector_processed


__all__ = [
    "create_land_ocean_mask_from_basins",
    "create_land_ocean_mask_from_shoreline",
    "remove_internal_polygon",
]
