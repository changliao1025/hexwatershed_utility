import os, sys

from pyearth.system.define_global_variables import *
from pyearth.toolbox.conversion.convert_vector_to_geojson import convert_vector_to_geojson
from pyearth.toolbox.conversion.convert_vector_to_global_raster import convert_vector_to_global_raster
from pyearth.toolbox.data.geoparquet.convert_geojson_to_geoparquet import convert_geojson_to_geoparquet

from pyearth.toolbox.management.vector.merge_files import merge_files
aRegion = list()
aRegion.append('af')
aRegion.append('ar')
aRegion.append('as')
aRegion.append('au')
aRegion.append('eu')
aRegion.append('gr')
aRegion.append('na')
aRegion.append('sa')
aRegion.append('si')

sWorkspace_data = '/compyfs/liao313/00raw/hydrology/hydroshed/hydrobasin'
sWorkspace_data2 = '/compyfs/liao313/00raw/hydrology/hydroshed/hydrobasin/geojson'
if not os.path.exists(sWorkspace_data2):
    os.makedirs(sWorkspace_data2)

aWatershed_boundary = list()
for sRegion in aRegion:
    sFolder = 'hybas_lake_' + sRegion + '_lev01-12_v1c'
    sFilename = 'hybas_lake_' + sRegion + '_lev03_v1c.shp'
    sWorkspace_hydroshed = sWorkspace_data + slash + sFolder
    sFilename_shp_in = os.path.join(sWorkspace_hydroshed, sFilename)
    sFilename_geojson_out = os.path.join(sWorkspace_data2, sFilename.replace('.shp', '.geojson'))
    #convert from shp to geojson
    #convert_vector_to_geojson(sFilename_shp_in, sFilename_geojson_out)
    aWatershed_boundary.append(sFilename_shp_in)


# Merge all the geojson files into a single file
sFilename_geojson_merge = sWorkspace_data2 + slash + 'hybas_lake_all_lev03.geojson'
if os.path.exists(sFilename_geojson_merge):
    os.remove(sFilename_geojson_merge)
    pass


merge_files(aWatershed_boundary, sFilename_geojson_merge)


convert_geojson_to_geoparquet(sFilename_geojson_merge, sFilename_geojson_merge.replace('.geojson', '.parquet'))

print('done')






