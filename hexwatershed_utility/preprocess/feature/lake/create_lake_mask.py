import os, sys

from pyearth.system.define_global_variables import *
from pyearth.toolbox.conversion.convert_vector_to_global_raster import convert_vector_to_global_raster



sFilename_in = '/compyfs/liao313/00raw/hydrology/hydroshed/hydrolake/global_lakes.geojson'
sFilename_tif_out = '/qfs/people/liao313/data/hexwatershed/global/raster/lakes_mask.tif'

#define resolution as 1km as the equator, which is
dResolution_x_in = 30 /3600.0
dResolution_y_in = 30 /3600.0


#covnert to raster

convert_vector_to_global_raster(sFilename_in, sFilename_tif_out, dResolution_x_in, dResolution_y_in,
                                iFlag_boundary_only_in = 0)

print('finished creating greatlakes mask in raster.')