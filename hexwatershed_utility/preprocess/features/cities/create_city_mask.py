import os, sys

from pyearth.system.define_global_variables import *


from pyearth.toolbox.conversion.convert_vector_to_global_raster import convert_vector_to_global_raster
sWorkspace_data = '/compyfs/liao313/00raw/city/simplemaps_worldcities_basicv1.77'
sFilename_geojson_out = os.path.join(sWorkspace_data, 'large_cities.geojson')

sFilename_tif_out = os.path.join(sWorkspace_data, 'large_cities.tif')
dResolution_x_in = 30 /3600.0
dResolution_y_in = 30 /3600.0

convert_vector_to_global_raster(sFilename_geojson_out, sFilename_tif_out, dResolution_x_in, dResolution_y_in)

print('finished creating watershed mask in raster.')