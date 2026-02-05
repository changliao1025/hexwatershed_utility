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

import os
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
from pyearthviz3d.geovista.utility import VisualizationConfig, AnimationConfig



def visualize_mpas_mesh(sFilename_mpas_mesh_in: str,
                          sFilename_out: Optional[str] = None,
                          dLongitude_focus_in: Optional[float] = 0.0,
                          dLatitude_focus_in: Optional[float] = 0.0,
                          dZoom_factor: float = 0.7,
                          iFlag_show_coastlines: bool = True,
                          iFlag_show_graticule: bool = True,
                          sCoastline_color: str = 'black',
                          dCoastline_width: float = 1.0,
                          iFlag_create_animation: Optional[bool] = False,
                         iAnimation_frames: Optional[int] = 36,
                         dAnimation_speed: Optional[float] = 1.0,
                         sAnimation_format: Optional[str] = 'mp4',
                         iFlag_wireframe_only: Optional[bool] = True,
                          dEdge_width: float = 1.0,
                          sEdge_color: str = 'black',
                         iFlag_cull_backfaces: Optional[bool] = True,
                         sCulling_mode: Optional[str] = 'auto',
                         iFlag_verbose_in: Optional[bool] = False) -> bool:

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
    # read the mpas mesh file
    import netCDF4 as nc
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

        dLongitude_focus_in = 0.0
        dLatitude_focus_in = 0.0
        sCoastline_color = "black"
        dCoastline_width = 1.0
        iFlag_verbose_in = 1

        pConfig = VisualizationConfig(
            longitude_focus=dLongitude_focus_in,
            latitude_focus=dLatitude_focus_in,
            window_size=(5000, 4000),
            show_coastlines=True,
            show_graticule=True,
            coastline_color=sCoastline_color,
            coastline_width=dCoastline_width,
            verbose=iFlag_verbose_in,
            )



        # Output or display

        map_single_frame(pMesh,
            aValid_cell_indices,
            pConfig,
            style = "wireframe",
            sScalar = None,
            sUnit = None,
            sFilename_out = sFilename_out,
            )

    except ImportError as e:
        logger.error('GeoVista library not available. Install with: pip install geovista')
        logger.error(f'Import error: {e}')
        return False

    except Exception as e:
        logger.error(f'Unexpected error during mesh visualization: {e}')
        logger.error(f'Error type: {type(e).__name__}')
        logger.error(f'Traceback: {traceback.format_exc()}')
        return False


