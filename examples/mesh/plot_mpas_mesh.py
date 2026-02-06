import os
from hexwatershed_utility.preprocess.mesh.visualize_mpas_mesh import visualize_mpas_mesh

import os, sys, platform

sPath_current = os.path.dirname(os.path.abspath(__file__))
sPath_library = os.path.dirname(os.path.dirname(sPath_current))
sys.path.append(sPath_library)




sFilename_mpas_mesh_in = os.path.join('/compyfs/liao313/04model/pyhexwatershed/global/pyflowline20260203001/jigsaw/out','culled_mesh.nc')
sFilename_mpas_mesh_in = 'C:\\scratch\\04model\\pyhexwatershed\\global\\pyflowline20260201001\\jigsaw\\out\\invert_mesh.nc'

sFilename_png_out = os.path.join('/qfs/people/liao313/workspace/python/hexwatershed_utility/figures', 'mpas_mesh.png')
sFilename_png_out = 'C:\\workspace\\python\\hexwatershed_utility\\figures\\mpas_mesh.jpg'
# Test the culling fix with explicit parameters
visualize_mpas_mesh(sFilename_mpas_mesh_in,
                   sFilename_out = None,
                   base_layer = 'natural_earth_1',
                   iFlag_wireframe_only = True,
                  window_size_in=(8000, 6000),  # Set the window size for better visibility
                   iFlag_verbose_in=True)  # Enable verbose output to see what's happening