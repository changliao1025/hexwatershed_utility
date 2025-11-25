
import os, sys, platform
from pyearth.visual.animate.animate_polyline_file_on_sphere import animate_polyline_file_on_sphere
from pyearth.visual.map.vector.map_vector_polyline_file import map_vector_polyline_file
sPlatform_os = platform.system()

if sPlatform_os == 'Windows':
    sPath  = 'C:\\workspace\\python\\hexwatershed_utility\\hexwatershed_utility'
    sys.path.append(os.path.dirname(sPath))

else:
    #macOS
    if sPlatform_os == 'Darwin':
        sFilename_source_mesh = '/Users/liao313/scratch/04model/pyhexwatershed/global/pyflowline20250927006//mpas.geojson' #use the L10-100 test mesh
        sFilename_polyline_in = '/Users/liao313/scratch/04model/pyhexwatershed/global/pyhexwatershed20250928001/hexwatershed/mpas_flow_direction.geojson'
    else:
        #linux
        pass

map_vector_polyline_file(sFilename_polyline_in,  iFlag_global_in=1,
                         sFilename_output_in='mpas_flow_direction_map.png')