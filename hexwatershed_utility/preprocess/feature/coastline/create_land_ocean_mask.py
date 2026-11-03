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


def create_land_ocean_mask_from_hydrobasins(
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


def simplify_vector_geometry(
    sFilename_vector_in: str,
    sFilename_vector_out: str,
    dTolerance: float = 0.001,
    verbose: bool = False,
) -> None:
    """
    Simplify polygon geometries in a vector file while preserving topology.

    Uses GDAL/OGR topology-preserving Douglas-Peucker simplification.
    """
    from osgeo import ogr
    from pyearth.gis.gdal.gdal_vector_format_support import get_vector_format_from_filename

    pInDs = ogr.Open(sFilename_vector_in, 0)
    if pInDs is None:
        raise RuntimeError(f"Could not open input vector dataset: {sFilename_vector_in}")

    pInLayer = pInDs.GetLayer()
    try:
        driver_name = get_vector_format_from_filename(sFilename_vector_out)
    except Exception:
        driver_name = "GeoJSON"
    pDriver = ogr.GetDriverByName(driver_name)
    if pDriver is None:
        pDriver = ogr.GetDriverByName("GeoJSON")

    if os.path.exists(sFilename_vector_out):
        try:
            os.remove(sFilename_vector_out)
        except Exception:
            pass

    pOutDs = pDriver.CreateDataSource(sFilename_vector_out)
    if pOutDs is None:
        pInDs = None
        raise RuntimeError(f"Failed to create output datasource: {sFilename_vector_out}")

    pOutLayer = pOutDs.CreateLayer(
        pInLayer.GetName(),
        srs=pInLayer.GetSpatialRef(),
        geom_type=pInLayer.GetGeomType(),
    )
    in_defn = pInLayer.GetLayerDefn()
    for i in range(in_defn.GetFieldCount()):
        pOutLayer.CreateField(in_defn.GetFieldDefn(i))

    pFeat = pInLayer.GetNextFeature()
    while pFeat:
        geom = pFeat.GetGeometryRef()
        if geom is not None and not geom.IsEmpty():
            simplified_geom = geom.SimplifyPreserveTopology(dTolerance)
            out_feat = ogr.Feature(pOutLayer.GetLayerDefn())
            out_feat.SetGeometry(simplified_geom if simplified_geom is not None else geom)
            for i in range(in_defn.GetFieldCount()):
                out_feat.SetField(i, pFeat.GetField(i))
            pOutLayer.CreateFeature(out_feat)
            out_feat = None
        pFeat = pInLayer.GetNextFeature()

    pInDs = None
    pOutDs = None


def create_land_ocean_mask_from_naturalearth(
    sWorkspace_output: str,
    sFilename_land_vector_in: Optional[str] = None,
    sResolution_nature_earth: str = "10m",
    iFlag_antarctic: int = 0,
    dResolution_x: float = 30 / 3600.0,
    dResolution_y: float = 30 / 3600.0,
    dThreshold_area_island: float = 1000.0,
    iRaster_buffer_pixel: int = 1,
    iFlag_remove_internal_rings: bool = False,
    iFlag_simplify_coastline: bool = False,
    dTolerance_simplify: float = 0.001,
    dDistance_buffer_meter: Optional[float] = None,
    iFlag_simplify_buffer: int = 1,
    sFilename_vector_out: Optional[str] = None,
    iFlag_convert_raster: bool = True,
    verbose: bool = False,
) -> Tuple[Optional[str], str]:
    """
    Create land-ocean mask from Natural Earth dataset or vector shoreline.

    Filters small islands, optionally removes internal lakes/holes, optionally simplifies
    coastlines, optionally creates geodesic buffer zones, and optionally rasterizes to GeoTIFF.

    Parameters
    ----------
    sWorkspace_output : str
        Directory to save outputs.
    sFilename_land_vector_in : Optional[str], default=None
        Path to land shoreline vector file. If None, auto-downloads Natural Earth land.
    sResolution_nature_earth : str, default='10m'
        Natural Earth resolution ('10m', '50m', or '110m') if downloading.
    iFlag_antarctic : int, default=0
        Include Antarctic coastline (0: exclude, 1: include).
    dResolution_x : float, default=30/3600.0
        X-resolution in degrees for rasterization.
    dResolution_y : float, default=30/3600.0
        Y-resolution in degrees for rasterization.
    dThreshold_area_island : float, default=1000.0
        Minimum area threshold in m² for island filtering.
    iRaster_buffer_pixel : int, default=1
        Number of pixels to buffer the raster coastline.
    iFlag_remove_internal_rings : bool, default=False
        If True, removes interior rings (inland lakes) from land polygons.
    iFlag_simplify_coastline : bool, default=False
        If True, simplifies coastline geometries using topology-preserving Douglas-Peucker.
    dTolerance_simplify : float, default=0.001
        Simplification tolerance in degrees.
    dDistance_buffer_meter : Optional[float], default=None
        Buffer distance in meters for geodesic vector buffering via pyearthbuffer.
    iFlag_simplify_buffer : int, default=1
        If 1, dissolve overlapping vector buffers.
    sFilename_vector_out : Optional[str], default=None
        Destination file path for output vector file. If None, auto-generated in sWorkspace_output.
    iFlag_convert_raster : bool, default=True
        If True, converts to global raster GeoTIFF. If False, skips rasterization.
    verbose : bool, default=False
        Verbosity level.

    Returns
    -------
    Tuple[Optional[str], str]
        (Path to output raster GeoTIFF or None, Path to cleaned/buffered vector file).
    """
    os.makedirs(sWorkspace_output, exist_ok=True)

    if sFilename_land_vector_in is None:
        sFilename_raw = os.path.join(sWorkspace_output, "land_ocean_mask_naturalearth.geojson")
        create_land_ocean_vector_mask_naturalearth(
            sFilename_raw,
            sResolution_coastal=sResolution_nature_earth,
            iFlag_exclude_antarctica=not bool(iFlag_antarctic),
            dThreshold_area_island=dThreshold_area_island,
            iFlag_remove_internal_rings=iFlag_remove_internal_rings,
        )
        sFilename_vector_cleaned = sFilename_raw
    else:
        sFilename_vector_cleaned = sFilename_land_vector_in
        if dThreshold_area_island is not None and dThreshold_area_island > 0:
            sFilename_wo_island = os.path.join(
                sWorkspace_output, "land_ocean_mask_wo_island.geojson"
            )
            remove_small_polygon(
                sFilename_vector_cleaned, sFilename_wo_island, dThreshold_area_island
            )
            sFilename_vector_cleaned = sFilename_wo_island

        if iFlag_remove_internal_rings:
            sFilename_wo_rings = os.path.join(
                sWorkspace_output, "land_ocean_mask_wo_internal_ring.geojson"
            )
            remove_internal_polygon(sFilename_vector_cleaned, sFilename_wo_rings)
            sFilename_vector_cleaned = sFilename_wo_rings

    # Coastline simplification
    if iFlag_simplify_coastline:
        sFilename_simplified = os.path.join(sWorkspace_output, "land_ocean_mask_simplified.geojson")
        simplify_vector_geometry(
            sFilename_vector_cleaned,
            sFilename_simplified,
            dTolerance=dTolerance_simplify,
            verbose=verbose,
        )
        sFilename_vector_current = sFilename_simplified
    else:
        sFilename_vector_current = sFilename_vector_cleaned

    # Vector buffering
    if dDistance_buffer_meter is not None and dDistance_buffer_meter > 0:
        if sFilename_vector_out is None:
            sFilename_vector_final = os.path.join(
                sWorkspace_output, "land_ocean_mask_buffered.geojson"
            )
        else:
            sFilename_vector_final = sFilename_vector_out

        gcsbuffer_file(
            input_file=sFilename_vector_current,
            output_file=sFilename_vector_final,
            distance=dDistance_buffer_meter,
            dissolve=bool(iFlag_simplify_buffer),
            min_area=dThreshold_area_island,
            verbose=verbose,
        )
    else:
        if sFilename_vector_out is not None and sFilename_vector_out != sFilename_vector_current:
            # Copy or convert to requested output format/destination
            from osgeo import ogr
            from pyearth.gis.gdal.gdal_vector_format_support import get_vector_format_from_filename

            pInDs = ogr.Open(sFilename_vector_current, 0)
            driver_name = get_vector_format_from_filename(sFilename_vector_out)
            pDriver = ogr.GetDriverByName(driver_name)
            if pDriver is not None and pInDs is not None:
                if os.path.exists(sFilename_vector_out):
                    try:
                        os.remove(sFilename_vector_out)
                    except Exception:
                        pass
                pDriver.CopyDataSource(pInDs, sFilename_vector_out)
                pInDs = None
            sFilename_vector_final = sFilename_vector_out
        else:
            sFilename_vector_final = sFilename_vector_current

    # Raster conversion
    sFilename_fixed: Optional[str] = None
    if iFlag_convert_raster:
        sFilename_tif_wo_island = os.path.join(sWorkspace_output, "land_ocean_mask_wo_island.tif")
        convert_vector_to_global_raster(
            sFilename_vector_current,
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

    return sFilename_fixed, sFilename_vector_final


__all__ = [
    "create_land_ocean_mask_from_hydrobasins",
    "create_land_ocean_mask_from_naturalearth",
    "simplify_vector_geometry",
]
