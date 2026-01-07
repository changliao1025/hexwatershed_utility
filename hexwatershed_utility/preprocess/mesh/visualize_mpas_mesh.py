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


    if not _validate_output_path(sFilename_out):
        return False

    # Create configuration object
    config = VisualizationConfig(
        longitude_focus=dLongitude_focus_in,
        latitude_focus=dLatitude_focus_in,
        zoom_factor=dZoom_factor,
        show_coastlines=iFlag_show_coastlines,
        show_graticule=iFlag_show_graticule,
        coastline_color=sCoastline_color,
        coastline_width=dCoastline_width,
        verbose=iFlag_verbose_in
    )

    animation_config = AnimationConfig(
        frames=iAnimation_frames,
        speed=dAnimation_speed,
        format=sAnimation_format
    ) if iFlag_create_animation else None

    try:
        # Import and setup GeoVista
        import geovista as gv
        if config.verbose:
            logger.info('Creating mesh visualization...')
            logger.info(f'  - Vertices: {len( aVertex_longititude)}')
            logger.info(f'  - Connectivity shape: { aConnectivity.shape}')
            logger.info(f'  - Focus: ({config.longitude_focus:.2f}°, {config.latitude_focus:.2f}°)')
            logger.info(f'  - Zoom factor: {config.zoom_factor}')

        # Validate connectivity array structure
        if aConnectivity.ndim != 2:
            logger.error(f'Connectivity array must be 2D, got {aConnectivity.ndim}D')
            return False

        if config.verbose:
            logger.info(f'Processing connectivity array of shape: {aConnectivity.shape}')
            logger.info(f'Number of vertices: {len(aVertex_longititude)}')

        # Convert 1-based connectivity indices to 0-based (MPAS uses 1-based indexing)
        # Invalid vertices are marked with fill_value (usually 0)
        fill_value = 0
        aConnectivity_0based = np.where(aConnectivity > fill_value, aConnectivity - 1, -1)

        # Count valid vertices per cell for diagnostics
        valid_vertices_per_cell = np.sum(aConnectivity_0based >= 0, axis=1)
        if config.verbose:
            logger.info(f'Valid vertices per cell - min: {np.min(valid_vertices_per_cell)}, '
                       f'max: {np.max(valid_vertices_per_cell)}, '
                       f'mean: {np.mean(valid_vertices_per_cell):.1f}')

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

        # Check for cells with too few vertices
        cells_with_few_vertices = np.sum(valid_vertices_per_cell < 3)
        if cells_with_few_vertices > 0:
            logger.warning(f'Found {cells_with_few_vertices} cells with fewer than 3 vertices')
            if config.verbose:
                # Log some examples
                bad_cells = np.where(valid_vertices_per_cell < 3)[0]
                for i in range(min(5, len(bad_cells))):
                    cell_idx = bad_cells[i]
                    logger.warning(f'  Cell {cell_idx}: {valid_vertices_per_cell[cell_idx]} vertices')

        # Transform to GeoVista unstructured mesh
        mesh = gv.Transform.from_unstructured(
             aVertex_longititude,
             aVertex_latitude,
            connectivity=connectivity_masked,
            crs=CRS
        )

        # Validate cell data array length matches mesh cells
        if len(aCellID) != mesh.n_cells:
            logger.error(f'Cell ID array length ({len( aCellID)}) does not match '
                        f'mesh cells ({mesh.n_cells})')
            return False

        # Prepare mesh metadata
        name = 'Mesh Cell ID'
        mesh.cell_data[name] =  aCellID

        # Setup plotter
        pPlotter = _setup_geovista_plotter(iFlag_off_screen=(sFilename_out is not None), iFlag_verbose_in=config.verbose)
        if pPlotter is None:
            return False

        if config.verbose:
            logger.info(f'Created GeoVista mesh with {mesh.n_cells} cells and {mesh.n_points} points')

        if iFlag_wireframe_only:
            # Wireframe mode: show only edges without filling (no colorbar needed)
            if config.verbose:
                logger.info('Rendering mesh in wireframe mode (edges only)')
                if iFlag_cull_backfaces:
                    logger.info(f'Face culling enabled with mode: {sCulling_mode}')

            # Determine culling mode
            if not iFlag_cull_backfaces:
                culling_mode = None
            elif sCulling_mode == 'auto':
                # For wireframe mode, front culling often works better
                culling_mode = 'front'
            elif sCulling_mode == 'none':
                culling_mode = None
            else:
                culling_mode = sCulling_mode

            # For wireframe with culling, we need a different approach
            # The issue is that transparent faces with culling can hide everything if normals are wrong
            if culling_mode is not None:
                try:
                    # Compute normals and ensure they point outward for sphere meshes
                    mesh = mesh.compute_normals(inplace=True, consistent=True, auto_orient_normals=True)
                    if config.verbose:
                        logger.info('Computed consistent face normals for proper culling')

                    # For wireframe mode, we'll use a hybrid approach:
                    # Show edges with very low opacity faces so culling works but faces are barely visible
                    pPlotter.add_mesh(
                        mesh,
                        show_edges=True,              # Show mesh edges
                        edge_color=sEdge_color,       # Color for mesh edges
                        line_width=dEdge_width,       # Width of edge lines
                        opacity=0.05,                 # Very low opacity so faces are barely visible but culling works
                        color='lightgray',            # Light color for barely visible faces
                        show_scalar_bar=False,        # No colorbar for wireframe-like mode
                        culling=culling_mode          # This works with low-opacity faces
                    )

                except Exception as e:
                    logger.warning(f'Could not compute face normals: {e}. Falling back to no culling.')
                    # Fall back to pure wireframe without culling
                    pPlotter.add_mesh(
                        mesh,
                        style='wireframe',
                        line_width=dEdge_width,
                        color=sEdge_color,
                        show_edges=True,
                        show_scalar_bar=False
                    )
            else:
                # No culling requested - use pure wireframe
                if config.verbose:
                    logger.info('Using pure wireframe mode (no culling)')
                pPlotter.add_mesh(
                    mesh,
                    style='wireframe',
                    line_width=dEdge_width,
                    color=sEdge_color,
                    show_edges=True,
                    show_scalar_bar=False
                )

            if config.verbose and culling_mode:
                logger.info(f'Using {culling_mode} culling for wireframe rendering with low-opacity faces')

        else:
            # Standard mode: show filled cells with scalars and colorbar
            if config.verbose and iFlag_cull_backfaces:
                logger.info(f'Face culling enabled with mode: {sCulling_mode}')

            sargs = {
                "title": name,
                "shadow": True,
                "title_font_size": 10,
                "label_font_size": 10,
                "fmt": "%.0f",  # Integer formatting for cell IDs
                "n_labels": 5,
            }

            # Determine culling mode
            if not iFlag_cull_backfaces:
                culling_mode = None
            elif sCulling_mode == 'auto':
                # For solid mode, back culling typically works best
                culling_mode = 'back'
            elif sCulling_mode == 'none':
                culling_mode = None
            else:
                culling_mode = sCulling_mode

            # Ensure proper face normals before culling
            if culling_mode is not None:
                try:
                    # Compute normals using PyVista's method
                    mesh = mesh.compute_normals(inplace=True)
                    if config.verbose:
                        logger.info('Computed face normals for proper culling')
                except Exception as e:
                    logger.warning(f'Could not compute face normals: {e}. Culling may not work properly.')
                    culling_mode = None

            if config.verbose and culling_mode:
                logger.info(f'Using {culling_mode} culling for solid rendering')

            pPlotter.add_mesh(
                mesh,
                scalars=name,
                scalar_bar_args=sargs,
                culling=culling_mode
            )

        # Configure camera
        _configure_camera(pPlotter, config)

        # Add geographic context
        _add_geographic_context(pPlotter, config)

        # Output or display
        return _handle_visualization_output(pPlotter, sFilename_out, config.verbose)

    except ImportError as e:
        logger.error('GeoVista library not available. Install with: pip install geovista')
        logger.error(f'Import error: {e}')
        return False

    except Exception as e:
        logger.error(f'Unexpected error during mesh visualization: {e}')
        logger.error(f'Error type: {type(e).__name__}')
        logger.error(f'Traceback: {traceback.format_exc()}')
        return False

def _handle_visualization_output(pPlotter, sFilename: Optional[str], iFlag_verbose_in: bool = False) -> bool:
    """
    Handle visualization output (save file or show interactive).

    Args:
        pPlotter: GeoVista plotter instance
        sFilename: Output filename or None for interactive
        iFlag_verbose_in: Enable verbose logging

    Returns:
        bool: True if successful, False otherwise
    """
    try:
        if sFilename is not None:
            # Save screenshot
            pPlotter.screenshot(sFilename)
            if iFlag_verbose_in:
                logger.info(f'✓ Visualization saved to: {sFilename}')

                # Verify file was created
                if os.path.exists(sFilename):
                    iFile_size = os.path.getsize(sFilename)
                    logger.info(f'  File size: {iFile_size / 1024:.1f} KB')
                else:
                    logger.warning(f'Screenshot command executed but file not found: {sFilename}')

            pPlotter.close()
            return True
        else:
            # Interactive display
            if iFlag_verbose_in:
                logger.info('Opening interactive visualization window...')
            pPlotter.show()
            return True

    except Exception as e:
        logger.error(f'Failed to handle visualization output: {e}')
        logger.error(f'Traceback: {traceback.format_exc()}')
        try:
            pPlotter.close()
        except Exception:
            pass
        return False



def _handle_single_frame_visualization(pMesh, sScalars: str, aValid_cell_indices: np.ndarray,
                                     sUnit: str, pConfig: VisualizationConfig,
                                     sFilename: Optional[str]) -> bool:
    """
    Handle single frame visualization (static image or interactive).

    Args:
        pMesh: GeoVista mesh object
        sScalars: Name of scalar field to visualize
        aValid_cell_indices: Indices of valid cells to display
        sUnit: Unit string for colorbar
        pConfig: Visualization configuration
        sFilename: Output filename or None for interactive

    Returns:
        bool: True if successful, False otherwise
    """
    try:
        # Setup plotter
        pPlotter = _setup_geovista_plotter(iFlag_off_screen=(sFilename is not None), iFlag_verbose_in=pConfig.verbose)
        if pPlotter is None:
            return False

        # Configure scalar bar
        dSargs = {
            "title": f"{sScalars} / {sUnit}" if sUnit else sScalars,
            "shadow": True,
            "title_font_size": 12,
            "label_font_size": 10,
            "fmt": "%.2f",
            "n_labels": 5,
        }

        # Add mesh to plotter - edges only with uniform color
        pMesh_valid = pMesh.extract_cells(aValid_cell_indices)

        # Ensure proper face normals for consistent rendering
        try:
            pMesh_valid = pMesh_valid.compute_normals(inplace=True)
        except Exception as e:
            logger.warning(f'Could not compute face normals: {e}')

        pPlotter.add_mesh(
            pMesh_valid,
            show_edges=True,          # Show mesh edges
            edge_color='black',         # Color for mesh edges
            line_width=1.0,           # Width of edge lines
            opacity=0.0,              # Make faces completely transparent
            show_scalar_bar=False,    # Hide colorbar
            culling='front'           # Use front culling for wireframe-like transparent faces
        )

        # Configure camera and add geographic context
        _configure_camera(pPlotter, pConfig)
        _add_geographic_context(pPlotter, pConfig)

        # Handle output
        return _handle_visualization_output(pPlotter, sFilename, pConfig.verbose)

    except Exception as e:
        logger.error(f'Error in single frame visualization: {e}')
        logger.error(f'Traceback: {traceback.format_exc()}')
        return False

def _handle_animation_visualization(pMesh, sScalars: str, aValid_cell_indices: np.ndarray,
                                  sUnit: str, pConfig: VisualizationConfig,
                                  pAnimation_config: AnimationConfig, sFilename: str) -> bool:
    """
    Handle animation visualization.

    Args:
        pMesh: GeoVista mesh object
        sScalars: Name of scalar field to visualize
        aValid_cell_indices: Indices of valid cells to display
        sUnit: Unit string for colorbar
        pConfig: Visualization configuration
        pAnimation_config: Animation configuration
        sFilename: Output animation filename

    Returns:
        bool: True if successful, False otherwise
    """
    if sFilename is None:
        logger.error('Animation mode requires output filename')
        return False

    try:
        # Setup off-screen plotter for animation
        pPlotter = _setup_geovista_plotter(iFlag_off_screen=True, iFlag_verbose_in=pConfig.verbose)
        if pPlotter is None:
            return False

        # Configure scalar bar
        dSargs = {
            "title": f"{sScalars} / {sUnit}" if sUnit else sScalars,
            "shadow": True,
            "title_font_size": 12,
            "label_font_size": 10,
            "fmt": "%.2f",
            "n_labels": 5,
        }

        # Add mesh to plotter - edges only with uniform color
        pMesh_valid = pMesh.extract_cells(aValid_cell_indices)

        # Ensure proper face normals for consistent rendering
        try:
            pMesh_valid = pMesh_valid.compute_normals(inplace=True)
        except Exception as e:
            logger.warning(f'Could not compute face normals: {e}')

        pPlotter.add_mesh(
            pMesh_valid,
            show_edges=True,          # Show mesh edges
            edge_color='red',         # Color for mesh edges
            line_width=2.0,           # Width of edge lines
            opacity=0.0,              # Make faces completely transparent
            show_scalar_bar=False,    # Hide colorbar
            culling='front'           # Use front culling for wireframe-like transparent faces
        )

        # Configure initial camera position
        _configure_camera(pPlotter, pConfig)

        # Add geographic context
        _add_geographic_context(pPlotter, pConfig)

        # Reset the zoom factor to 1.0 so it won't zoom in too much during the animation
        pConfig.zoom_factor = 1.0

        # Create animation
        iFlag_success = _create_rotation_animation(
            pPlotter, sFilename, pConfig.longitude_focus, pConfig.latitude_focus,
            pConfig.zoom_factor, pAnimation_config.frames, pAnimation_config.speed,
            pAnimation_config.format, pConfig.verbose
        )

        pPlotter.close()
        return iFlag_success

    except Exception as e:
        logger.error(f'Error in animation visualization: {e}')
        logger.error(f'Traceback: {traceback.format_exc()}')
        return False

def _create_rotation_animation(pPlotter, sFilename_out, dLongitude_start, dLatitude_focus,
                               dZoom_factor, iAnimation_frames, dAnimation_speed, sAnimation_format, iFlag_verbose_in):
    """
    Create a rotating animation of the 3D globe visualization.

    Generates multiple frames by rotating the camera around the globe with enhanced
    camera movement patterns, then combines them into a video file.

    Args:
        pPlotter: GeoVista plotter instance with mesh already added
        sFilename_out (str): Output animation file path (e.g., 'animation.mp4')
        dLongitude_start (float): Starting longitude for rotation in degrees
        dLatitude_focus (float): Base latitude for camera focus in degrees
        dZoom_factor (float): Camera zoom level
        iAnimation_frames (int): Number of frames for 360° rotation
        dAnimation_speed (float): Degrees per frame
        sAnimation_format (str): Output format ('mp4', 'gif', 'avi')
        iFlag_verbose_in (bool): Enable verbose logging

    Returns:
        bool: True if animation created successfully, False otherwise
    """
    try:
        if iFlag_verbose_in:
            logger.info(f'Creating {iAnimation_frames}-frame rotation animation')
            logger.info(f'  - Starting longitude: {dLongitude_start:.1f}°')
            logger.info(f'  - Base latitude: {dLatitude_focus:.1f}°')
            logger.info(f'  - Rotation speed: {dAnimation_speed:.1f}°/frame')
            logger.info(f'  - Output format: {sAnimation_format}')

        # Animation parameters
        dEarth_radius = DEFAULT_EARTH_RADIUS
        dCamera_distance = dEarth_radius * DEFAULT_CAMERA_DISTANCE_MULTIPLIER
        dAmplitude_deg = 20.0  # Latitude oscillation amplitude
        dCycles = 1.0  # Number of sine cycles over full rotation
        dPhase = 0.0  # Phase shift for sine wave

        # Initialize movie recording
        pPlotter.open_movie(sFilename_out, framerate=30)

        if iFlag_verbose_in:
            logger.info('Generating animation frames...')

        for iFrame in range(iAnimation_frames):
            # Calculate current longitude with smooth rotation
            dLongitude_current = dLongitude_start + (iFrame * dAnimation_speed)
            dLongitude_current = dLongitude_current % 360.0  # Keep within [0, 360)
            if dLongitude_current > 180.0:
                dLongitude_current -= 360.0  # Convert to [-180, 180]

            # Enhanced latitude movement: sine-wave oscillation for dynamic viewing
            # This creates a more interesting camera path than fixed latitude
            dFrames_div = float(iAnimation_frames) if iAnimation_frames > 0 else 1.0
            dTheta = 2.0 * math.pi * (float(iFrame) / dFrames_div) * dCycles + dPhase
            dLatitude_current = float(dLatitude_focus) + dAmplitude_deg * math.sin(dTheta)

            # Clamp latitude to avoid pole singularities
            dLatitude_current = max(-89.9, min(89.9, dLatitude_current))

            # Convert to radians for calculations
            dLon_rad = math.radians(dLongitude_current)
            dLat_rad = math.radians(dLatitude_current)

            # Calculate focal point on Earth surface
            dX_focal = dEarth_radius * math.cos(dLat_rad) * math.cos(dLon_rad)
            dY_focal = dEarth_radius * math.cos(dLat_rad) * math.sin(dLon_rad)
            dZ_focal = dEarth_radius * math.sin(dLat_rad)

            # Calculate camera position away from Earth
            dX_camera = dCamera_distance * math.cos(dLat_rad) * math.cos(dLon_rad)
            dY_camera = dCamera_distance * math.cos(dLat_rad) * math.sin(dLon_rad)
            dZ_camera = dCamera_distance * math.sin(dLat_rad)

            aFocal_point = [dX_focal, dY_focal, dZ_focal]
            aCamera_position = [dX_camera, dY_camera, dZ_camera]

            # Update camera with smooth transitions
            pPlotter.camera.focal_point = aFocal_point
            pPlotter.camera.position = aCamera_position
            pPlotter.camera.up = [0, 0, 1]  # Maintain Z-up orientation

            # Apply zoom factor for consistent view
            pPlotter.camera.zoom(dZoom_factor)

            # Ensure axes remain visible throughout animation
            try:
                pPlotter.add_axes()
            except Exception:
                pass  # Axes may already exist

            # Render the current frame
            pPlotter.render()

            try:
                pPlotter.write_frame()

                if iFlag_verbose_in and (iFrame + 1) % max(1, iAnimation_frames // 10) == 0:
                    dProgress = ((iFrame + 1) / iAnimation_frames) * 100
                    logger.info(f'  Progress: {dProgress:.0f}% ({iFrame + 1}/{iAnimation_frames} frames)')

            except Exception as e:
                logger.error(f'Failed to render frame {iFrame + 1}: {e}')
                try:
                    pPlotter.close()
                except Exception:
                    pass
                return False

        # Close movie recording
        try:
            pPlotter.close()
        except Exception as e:
            logger.warning(f'Error closing plotter: {e}')

        # Validate output file creation
        if not os.path.exists(sFilename_out):
            logger.error('Animation file was not created')
            return False

        # Log success information
        iFile_size = os.path.getsize(sFilename_out)
        if iFlag_verbose_in:
            logger.info(f'✓ Animation created successfully: {sFilename_out}')
            logger.info(f'  File size: {iFile_size / (1024*1024):.2f} MB')
            logger.info(f'  Frames: {iAnimation_frames}')
            logger.info(f'  Format: {sAnimation_format.upper()}')
            logger.info(f'  Duration: ~{iAnimation_frames / 30:.1f} seconds at 30 FPS')

        return True

    except Exception as e:
        logger.error(f'Unexpected error during animation creation: {e}')
        logger.error(f'Traceback: {traceback.format_exc()}')
        try:
            pPlotter.close()
        except Exception:
            pass
        return False

