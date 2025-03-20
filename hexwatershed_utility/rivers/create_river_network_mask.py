import os, sys

from pyearth.system.define_global_variables import *
from pyearth.toolbox.conversion.convert_vector_to_global_raster import convert_vector_to_global_raster

sFilename_vector_out= '/qfs/people/liao313/data/hexwatershed/global/vector/river_networks_4.0E4_wo_lakes.geojson'

sFilename_tif_out = '/qfs/people/liao313/data/hexwatershed/global/raster/global_river_networks_4.0E4.tif'

#define resolution as 1km as the equator, which is
dResolution_x_in = 30 /3600.0
dResolution_y_in = 30 /3600.0


#covnert to raster

convert_vector_to_global_raster(sFilename_vector_out, sFilename_tif_out, dResolution_x_in, dResolution_y_in)

print('finished creating river network raster.')