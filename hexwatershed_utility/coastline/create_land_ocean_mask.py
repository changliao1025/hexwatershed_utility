import os, sys
from pyearth.system.define_global_variables import *

from pyearth.toolbox.data.ocean.define_land_ocean_mask import create_land_ocean_vector_mask

from pyearth.toolbox.conversion.convert_vector_to_global_raster import convert_vector_to_global_raster
from pyearth.toolbox.management.vector.remove_small_polygon import remove_small_polygon

sFilename_geojson_raw = '/qfs/people/liao313/data/hexwatershed/global/vector/land_ocean_mask.geojson'
sWorkspace_vector = '/qfs/people/liao313/data/hexwatershed/global/vector/'
sWorkspace_raster = '/qfs/people/liao313/data/hexwatershed/global/raster/'

#define resolution as 1km as the equator, which is
dResolution_x_in = 30 /3600.0
dResolution_y_in = 30 /3600.0

#create_land_ocean_vector_mask(sFilename_geojson_out )

iFlag_remove_small_island = 1
iFlag_simplify_coastline = 1

dThroshold_area = 3.0E4  #in square km for small islans
sTolerance_island = '3E4'
dSimplify_tolerance = 1.0E6 #in meters
sTolerance_distance = '1E6'

if iFlag_remove_small_island == 1:
    sFilename_geojson_in = sFilename_geojson_raw
    if iFlag_simplify_coastline == 1:
        sTemp = 'land_ocean_mask_wo_island_' + sTolerance_island +'_' +sTolerance_distance + '.geojson'
        sFilename_geojson_out = os.path.join(sWorkspace_vector, sTemp)
        remove_small_polygon(sFilename_geojson_in, sFilename_geojson_out, dThroshold_area,
                          dSimplify_tolerance = dSimplify_tolerance)
        sTemp = 'land_ocean_mask_wo_land_' + sTolerance_island +'_' + sTolerance_distance + '.tif'
        sFilename_tif_out = os.path.join(sWorkspace_raster, sTemp)
        print('converting to raster', sFilename_tif_out)
        convert_vector_to_global_raster(sFilename_geojson_out, sFilename_tif_out,
                                         dResolution_x_in, dResolution_y_in,
                                        iFlag_boundary_only_in = 0 )
    else:
        sTemp = 'land_ocean_mask_wo_island_' + sTolerance_island +'.geojson'
        sFilename_geojson_out = os.path.join(sWorkspace_vector, sTemp)
        remove_small_polygon(sFilename_geojson_in, sFilename_geojson_out, dThroshold_area   )
        sTemp = 'land_ocean_mask_wo_island_' + sTolerance_island +'.tif'
        sFilename_tif_out = os.path.join(sWorkspace_raster, sTemp)
        print('converting to raster', sFilename_tif_out)
        convert_vector_to_global_raster(sFilename_geojson_out, sFilename_tif_out,
                                         dResolution_x_in, dResolution_y_in,
                                        iFlag_boundary_only_in = 0 )
else:
    sFilename_geojson_out = sFilename_geojson_raw
    print('finished creating land ocean mask in vector.')
    #covnert to raster
    sTemp = 'land_ocean_mask.tif'
    sFilename_tif_out = os.path.join(sWorkspace_raster, sTemp)
    print('converting to raster', sFilename_tif_out)
    convert_vector_to_global_raster(sFilename_geojson_out, sFilename_tif_out,
                                     dResolution_x_in, dResolution_y_in,
                                    iFlag_boundary_only_in = 0 )

print('Finished creating land ocean mask in raster.')