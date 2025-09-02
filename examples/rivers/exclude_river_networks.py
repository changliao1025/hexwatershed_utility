import os, sys

from pyearth.system.define_global_variables import *
from pyearth.toolbox.conversion.convert_vector_to_global_raster import convert_vector_to_global_raster


from pyearth.toolbox.analysis.extract.exclude_vector_by_polygon_files import exclude_vector_by_polygon_files


sFilename_river_network_out = '/qfs/people/liao313/data/hexwatershed/conus/vector/river_networks_wo_greatlakes.geojson'

sFilename_vector_in = '/compyfs/liao313/04model/pyflowline/global/flowline_hydroshed_simplified_1.0E4.geojson'
sFilename_vector_in = '/compyfs/liao313/04model/pyflowline/global/flowline_hydroshed_simplified_2.0E4.geojson'
sFilename_vector_in = '/compyfs/liao313/04model/pyflowline/global/flowline_hydroshed_simplified_4.0E4.geojson'

sWorkspace_data_out  = '/compyfs/liao313/00raw/hydrology/hydroshed/hydrolake/'


sFilename_out = '/compyfs/liao313/00raw/hydrology/hydroshed/hydrolake/global_lakes.geojson'
aFilename_polygon_in=[sFilename_out]
sFilename_vector_out= '/qfs/people/liao313/data/hexwatershed/global/vector/river_networks_4.0E4_wo_lakes.geojson'
exclude_vector_by_polygon_files(sFilename_vector_in, aFilename_polygon_in, sFilename_vector_out)
print('finished excluding lakes from river network.')
