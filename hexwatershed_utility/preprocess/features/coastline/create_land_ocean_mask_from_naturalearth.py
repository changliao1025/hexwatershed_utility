import os
import glob
from osgeo import gdal, ogr
gdal.UseExceptions()

from pyearth.toolbox.management.vector.remove_small_polygon import remove_small_polygon
from pyearth.toolbox.management.vector.merge_features import merge_features
from pyearth.toolbox.conversion.convert_vector_to_global_raster import convert_vector_to_global_raster
from pyearth.toolbox.data.ocean.define_land_ocean_mask import create_land_ocean_vector_mask_naturalearth
from pyearth.toolbox.analysis.image.raster_process import fix_raster_antimeridian_issue
from pyearth.toolbox.analysis.image.raster_process import create_raster_buffer_zone, fix_raster_antimeridian_issue
from pyearthbuffer.utility.create_gcs_buffer_zone import create_buffer_zone_polygon_file
def create_land_ocean_mask_from_naturalearth(sWorkspace_coastline_output,
                                                                             dResolution_x_in, dResolution_y_in,
                                                                             dThreshold_area_island,
                                                                             dResolution_coastline_buffer,
                                                                             iRaster_buffer_pixel = 1):

    sFilename_naturalearth = os.path.join(sWorkspace_coastline_output, 'land_ocean_mask_naturalearth.geojson')
    create_land_ocean_vector_mask_naturalearth(sFilename_naturalearth)

    sFilename_wo_island = os.path.join(sWorkspace_coastline_output, 'land_ocean_mask_wo_island.geojson')
    remove_small_polygon(sFilename_naturalearth, sFilename_wo_island, dThreshold_area_island )

    sFilename_tif_wo_island = os.path.join(sWorkspace_coastline_output, 'land_ocean_mask_wo_island.tif')
    convert_vector_to_global_raster(sFilename_wo_island,
                                    sFilename_tif_wo_island,
                                    dResolution_x_in,
                                    dResolution_y_in,
                                    iFlag_boundary_only_in = 0,
                                    dFill_value_in = 2)


    if iRaster_buffer_pixel > 0:
        sFilename_tif_wo_island_buffered = os.path.join(sWorkspace_coastline_output, 'land_ocean_mask_wo_island_buffered.tif')
        create_raster_buffer_zone(sFilename_tif_wo_island, sFilename_tif_wo_island_buffered, 1, int(iRaster_buffer_pixel) )
        sFilename_tif_wo_island_buffered_fixed = os.path.join(sWorkspace_coastline_output, 'land_ocean_mask_wo_island_buffered_fixed.tif')
        #fix dateline issue
        fix_raster_antimeridian_issue(sFilename_tif_wo_island_buffered, sFilename_tif_wo_island_buffered_fixed, 1,2, iRaster_buffer_pixel )
    else:
        sFilename_tif_wo_island_buffered_fixed = os.path.join(sWorkspace_coastline_output, 'land_ocean_mask_wo_island_buffered_fixed.tif')
        fix_raster_antimeridian_issue(sFilename_tif_wo_island_buffered, sFilename_tif_wo_island_buffered_fixed, 1,2, 1)

    return sFilename_tif_wo_island_buffered_fixed, sFilename_wo_island
