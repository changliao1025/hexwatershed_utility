import os
from hexwatershed_utility.preprocess.mesh.visualize_mpas_mesh import visualize_mpas_mesh

import os, sys, platform

sPath_current = os.path.dirname(os.path.abspath(__file__))
sPath_library = os.path.dirname(os.path.dirname(sPath_current))
sys.path.append(sPath_library)




sFilename_mpas_mesh_in = os.path.join('/compyfs/liao313/04model/pyhexwatershed/global/pyflowline20260203001/jigsaw/out','culled_mesh.nc')

sFilename_png_out = os.path.join('/qfs/people/liao313/workspace/python/hexwatershed_utility/figures', 'mpas_mesh.png')
# Test the culling fix with explicit parameters
visualize_mpas_mesh(sFilename_mpas_mesh_in,
                   sFilename_out = sFilename_png_out,
                   iFlag_wireframe_only=True,
                   iFlag_cull_backfaces=True,
                   sCulling_mode='back',  # Use auto mode for best results
                   iFlag_verbose_in=True)  # Enable verbose output to see what's happening