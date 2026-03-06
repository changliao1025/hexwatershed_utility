"""
Visualization module for uraster class.

This module contains all visualization-related methods that were moved from the main uraster class
to reduce the size of the main uraster.py file and improve code organization.

Features:
- 3D mesh visualization using GeoVista
- Interactive and static rendering modes
- Animation support with rotation and camera movement
- Comprehensive error handling and validation
- Support for multiple output formats
"""

import os,sys
import logging
import traceback
import math
from typing import Optional, List, Tuple, Union, Dict, Any
import numpy as np
from osgeo import gdal, ogr
gdal.UseExceptions()

# Set up logging
logger = logging.getLogger(__name__)
CRS = "EPSG:4326"

# Constants for visualization
DEFAULT_EARTH_RADIUS = 1.0
DEFAULT_CAMERA_DISTANCE_MULTIPLIER = 3.0
DEFAULT_ZOOM_FACTOR = 0.7
VALID_ANIMATION_FORMATS = ['mp4', 'gif', 'avi']
VALID_IMAGE_FORMATS = ['.png', '.jpg', '.jpeg', '.svg', '.tif', '.tiff']
COORDINATE_BOUNDS = {'longitude': (-180, 180), 'latitude': (-90, 90)}

from pyearthviz3d.geovista.map_single_frame import map_single_frame
from pyearthviz3d.geovista.animate_rotating_frames import animate_rotating_frames
from pyearthviz3d.geovista.utility import VisualizationConfig, AnimationConfig, ScalarBarConfig

def map_mpas_mesh_on_sphere(sFilename_mpas_mesh_in: str,
                          sFilename_out: Optional[str] = None,
    **kwargs,) -> bool:

    """
    Visualize the source mesh topology using GeoVista 3D globe rendering.

    Creates an interactive or saved 3D visualization of the unstructured mesh
    with proper geographic context including coastlines and coordinate grid.

    Args:
        sFilename_out: Output screenshot file path. If None, displays interactive viewer.
            Supports formats: .png, .jpg, .svg
        dLongitude_focus_in: Camera focal point longitude in degrees (-180 to 180).
            Default is 0.0 (prime meridian).
        dLatitude_focus_in: Camera focal point latitude in degrees (-90 to 90).
            Default is 0.0 (equator).
        dZoom_factor: Camera zoom level. Higher values zoom in. Default is 0.7.
        iFlag_show_coastlines: Show coastline overlay. Default is True.
        iFlag_show_graticule: Show coordinate grid with labels. Default is True.
        sCoastline_color: Color for coastlines. Default is 'black'.
            Examples: 'white', 'red', 'blue', 'gray', or RGB tuples like (1.0, 0.0, 0.0).
        dCoastline_width: Line width for coastlines. Default is 1.0.
        iFlag_cull_backfaces: Hide back faces of the sphere. Default is True.
            When True, only faces facing the camera are rendered, preventing see-through effect.
        sCulling_mode: Culling strategy ('auto', 'back', 'front', 'none'). Default is 'auto'.
            'auto' selects the best culling mode based on rendering style.
            'back' culls back-facing polygons (traditional for solid surfaces).
            'front' culls front-facing polygons (good for wireframes).
            'none' disables culling (shows all faces).
        iFlag_verbose_in: If True, print detailed progress messages. Default is False.

    Returns:
        True if visualization successful, False otherwise

    Note:
        - Requires 'geovista' package: pip install geovista
        - Interactive mode requires display environment
        - Mesh topology must be built before visualization (call rebuild_mesh_topology first)
    """
    defaults = {
        "dLongitude_focus_in": 0.0,
        "dLatitude_focus_in": 0.0,
        "dImage_scale_in": 1.0,
        "dZoom_factor": 0.7,
        "window_size_in": (800, 600),
        "iFlag_show_coastlines": True,
        "iFlag_show_graticule": True,
        "sColormap": "viridis",
        "sCoastline_color": "black",
        "sBase_layer": None, #'natural_earth_1',
        "dCoastline_width": 1.0,
        "iFlag_wireframe_only": False,
        "iFlag_create_animation": False,
        "iAnimation_frames": 36,
        "dAnimation_speed": 1.0,
        "sAnimation_format": "mp4",
        "dEdge_width": 1.0,
        "sEdge_color": "black",
        "iFlag_verbose_in": False,
        "sVariable_to_plot": None,
    }
    # Merge defaults with provided kwargs
    merged_params = {**defaults, **kwargs}

    # Extract parameters
    dLongitude_focus_in = merged_params["dLongitude_focus_in"]
    dLatitude_focus_in = merged_params["dLatitude_focus_in"]
    dImage_scale_in = merged_params["dImage_scale_in"]
    dZoom_factor = merged_params["dZoom_factor"]
    window_size_in = merged_params["window_size_in"]
    iFlag_show_coastlines = merged_params["iFlag_show_coastlines"]
    iFlag_show_graticule = merged_params["iFlag_show_graticule"]
    sColormap = merged_params["sColormap"]
    sCoastline_color = merged_params["sCoastline_color"]
    sBase_layer = merged_params["sBase_layer"]
    dCoastline_width = merged_params["dCoastline_width"]
    iFlag_create_animation = merged_params["iFlag_create_animation"]
    iFlag_wireframe_only = merged_params["iFlag_wireframe_only"]
    iAnimation_frames = merged_params["iAnimation_frames"]
    dAnimation_speed = merged_params["dAnimation_speed"]
    sAnimation_format = merged_params["sAnimation_format"]
    iFlag_verbose_in = merged_params["iFlag_verbose_in"]
    dEdge_width = merged_params["dEdge_width"]
    sEdge_color = merged_params["sEdge_color"]
    sVariable_to_plot = merged_params["sVariable_to_plot"]
    # read the mpas mesh file
    import netCDF4 as nc
    # Define a base window size and font size
    base_window_size = 1200
    base_font_size = 8
    base_line_width = 0.1

    # Calculate scaling factor
    scaling_factor = max(window_size_in) / base_window_size

    # Adjust font size
    scaled_font_size = int(base_font_size * scaling_factor)
    scaled_line_width = base_line_width * scaling_factor
    try:
        pDataset = nc.Dataset(sFilename_mpas_mesh_in, 'r')
    except Exception as e:
        logger.error(f'Could not open MPAS mesh file {sFilename_mpas_mesh_in}: {e}')
        return False
    try:
        #read mesh lon and lat, convert to degrees
        aVertex_longititude = pDataset.variables['lonVertex'][:] * (180.0 / np.pi)
        #fix range of longitudes if necessary from [0, 360] to [-180, 180]
        aVertex_longititude = np.where(aVertex_longititude > 180, aVertex_longititude - 360, aVertex_longititude)
        aVertex_latitude = pDataset.variables['latVertex'][:] * (180.0 / np.pi)
        aConnectivity = pDataset.variables['verticesOnCell'][:]
        aCellID = pDataset.variables['indexToCellID'][:]

    except Exception as e:
        logger.error(f'Error reading mesh variables from {sFilename_mpas_mesh_in}: {e}')
        pDataset.close()
        return False

    try:
        # Import and setup GeoVista
        import geovista as gv
        # Validate connectivity array structure
        if aConnectivity.ndim != 2:
            logger.error(f'Connectivity array must be 2D, got {aConnectivity.ndim}D')
            return False

        # Convert 1-based connectivity indices to 0-based (MPAS uses 1-based indexing)
        # Invalid vertices are marked with fill_value (usually 0)
        fill_value = 0
        aConnectivity_0based = np.where(aConnectivity > fill_value, aConnectivity - 1, -1)

        # Create masked connectivity array (mask invalid indices)
        connectivity_masked = np.ma.masked_where(
            aConnectivity_0based < 0,
            aConnectivity_0based
        )

        # Validate connectivity indices (after conversion to 0-based)
        valid_connectivity = aConnectivity_0based[aConnectivity_0based >= 0]
        if len(valid_connectivity) > 0:
            max_vertex_idx = len(aVertex_longititude) - 1
            if np.max(valid_connectivity) > max_vertex_idx:
                logger.error(f'Connectivity contains invalid vertex index: '
                           f'max={np.max(valid_connectivity)}, vertices={len(aVertex_longititude)}')
                return False

        # Transform to GeoVista unstructured mesh
        pMesh = gv.Transform.from_unstructured(
             aVertex_longititude,
             aVertex_latitude,
            connectivity=connectivity_masked,
            crs=CRS
        )

        aValid_cell_indices = np.arange(len(aCellID))

        if sVariable_to_plot is not None and sVariable_to_plot in pDataset.variables:
            try:
                variable_data = pDataset.variables[sVariable_to_plot][:]
                if variable_data.shape[0] != len(aCellID):
                    logger.error(f'Variable {sVariable_to_plot} length does not match number of cells')
                    return False
                pMesh.cell_data[sVariable_to_plot] = np.sqrt(variable_data)
                sScalar = sVariable_to_plot
            except Exception as e:
                logger.error(f'Error reading variable {sVariable_to_plot} from dataset: {e}')
                return False
        else:
            sScalar = "Cell ID"
            pMesh.cell_data[sScalar] = aCellID

        config_static = VisualizationConfig(
            longitude_focus=dLongitude_focus_in,
            latitude_focus=dLatitude_focus_in,
            image_scale=dImage_scale_in,
            window_size=window_size_in,
            zoom_factor=dZoom_factor,
            show_coastlines=iFlag_show_coastlines,
            show_graticule=iFlag_show_graticule,
            colormap=sColormap,
            coastline_color=sCoastline_color,
            coastline_width=dCoastline_width,
            verbose=iFlag_verbose_in,
            base_layer=sBase_layer,
        )
        config_anima = (
            AnimationConfig(
                frames=iAnimation_frames,
                speed=dAnimation_speed,
                format=sAnimation_format,
                longitude_start=dLongitude_focus_in,
                latitude_start=dLatitude_focus_in,
            )
            if iFlag_create_animation
            else None
        )
        scalar_config = ScalarBarConfig(
                    title=sScalar,
                    title_font_size=scaled_font_size,
                    label_font_size=scaled_font_size,
                    position_x=0.8,
                    position_y=0.25,
                    orientation="vertical",
                    n_labels=5,
                )

        config_colorbar = scalar_config if sVariable_to_plot is not None else None

        if iFlag_wireframe_only:
            style = "wireframe"
            if config_static.verbose:
                logger.info("Using wireframe-only visualization mode")
        else:
            style = "surface"
            if config_static.verbose:
                logger.info("Using surface visualization mode")

        # Handle animation vs single frame visualization
        if config_anima is not None:
            animate_rotating_frames(
                pMesh,
                aValid_cell_indices,
                config_static,
                config_anima,
                style = style,
                sFilename_out=sFilename_out,
                scalar_config=config_colorbar,
            )
        else:
            map_single_frame(
                pMesh,
                aValid_cell_indices,
                config_static,
                sScalar = sScalar,
                style = style,
                sFilename_out=sFilename_out,
                scalar_config=config_colorbar,
            )

        return True

    except ImportError as e:
        logger.error('GeoVista library not available. Install with: pip install geovista')
        logger.error(f'Import error: {e}')
        return False

    except Exception as e:
        logger.error(f'Unexpected error during mesh visualization: {e}')
        logger.error(f'Error type: {type(e).__name__}')
        logger.error(f'Traceback: {traceback.format_exc()}')
        return False


