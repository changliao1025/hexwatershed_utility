import os, sys
from pyearth.system.define_global_variables import *

from pyearth.toolbox.data.ocean.define_land_ocean_mask import create_land_ocean_vector_mask

from pyearth.toolbox.conversion.convert_vector_to_global_raster import convert_vector_to_global_raster
from pyearth.toolbox.management.vector.remove_small_polygon import remove_small_polygon
sPath_project = '/qfs/people/liao313/workspace/python/hexwatershed_utility/'
#add the project path of the pythonpath
sys.path.append(sPath_project)
from hexwatershed_utility.codes.merge_features import merge_features
sFilename_geojson_raw = '/compyfs/liao313/00raw/mesh/conus/ICoM_CONUS/domain_CONUS.geojson'
sWorkspace_vector = '/qfs/people/liao313/data/hexwatershed/global/vector/'
sWorkspace_raster = '/qfs/people/liao313/data/hexwatershed/global/raster/'

#define resolution as 1km as the equator, which is
dResolution_x_in = 30 /3600.0
dResolution_y_in = 30 /3600.0

dResolution_x_in = 1.0/8
dResolution_y_in = 1.0/8
sResolution = "{:.2E}".format(dResolution_x_in)

#create_land_ocean_vector_mask(sFilename_geojson_out )

iFlag_remove_small_island = 1
iFlag_simplify_coastline = 0

dThreshold_area = 3.0E4  #in square km for small islans
sThreshold_island = "{:.2E}".format(dThreshold_area)
dSimplify_tolerance = 1.0E6 #in meters
sSimplify_tolerance =  "{:.2E}".format(dSimplify_tolerance)

if iFlag_remove_small_island == 1:
    sFilename_geojson_in = sFilename_geojson_raw
    if iFlag_simplify_coastline == 1:
        sTemp = 'land_ocean_mask_wo_island_' + sThreshold_island +'_' +sSimplify_tolerance + '.geojson'
        sFilename_geojson_out = os.path.join(sWorkspace_vector, sTemp)
        remove_small_polygon(sFilename_geojson_in, sFilename_geojson_out, sThreshold_island,
                          dSimplify_tolerance = dSimplify_tolerance)
        sTemp = 'land_ocean_mask_wo_island_' + sThreshold_island +'_' + sSimplify_tolerance+'_' + sResolution + '.tif'
        sFilename_tif_out = os.path.join(sWorkspace_raster, sTemp)
        print('converting to raster', sFilename_tif_out)
        convert_vector_to_global_raster(sFilename_geojson_out, sFilename_tif_out,
                                         dResolution_x_in, dResolution_y_in,
                                        iFlag_boundary_only_in = 0 )
    else:
        sTemp = 'land_ocean_mask_wo_island_' + sThreshold_island + '_conus' + '.geojson'
        sFilename_geojson_out = os.path.join(sWorkspace_vector, sTemp)
        remove_small_polygon(sFilename_geojson_in, sFilename_geojson_out, sThreshold_island   )
        sTemp = 'land_ocean_mask_wo_island_' + sThreshold_island + '_conus_merged' + '.geojson'
        sFilename_clip_new = os.path.join(sWorkspace_vector, sTemp)
        merge_features(sFilename_geojson_out, sFilename_clip_new)
        sTemp = 'land_ocean_mask_wo_island_' + sThreshold_island +'_' + sResolution + '_conus' +'.tif'
        sFilename_tif_out = os.path.join(sWorkspace_raster, sTemp)
        print('converting to raster', sFilename_tif_out)
        convert_vector_to_global_raster(sFilename_clip_new, sFilename_tif_out,
                                         dResolution_x_in, dResolution_y_in,
                                        iFlag_boundary_only_in = 0 )
else:
    sFilename_geojson_out = sFilename_geojson_raw
    print('finished creating land ocean mask in vector.')
    #covnert to raster
    sTemp = 'land_ocean_mask' + '_' + sResolution + '.tif'
    sFilename_tif_out = os.path.join(sWorkspace_raster, sTemp)
    print('converting to raster', sFilename_tif_out)
    convert_vector_to_global_raster(sFilename_geojson_out, sFilename_tif_out,
                                     dResolution_x_in, dResolution_y_in,
                                    iFlag_boundary_only_in = 0 )

print('Finished creating land ocean mask in raster.')