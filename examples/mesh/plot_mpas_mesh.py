import os
from hexwatershed_utility.preprocess.mesh.visualize_mpas_mesh import visualize_mpas_mesh

import os, sys, platform

sPath_current = os.path.dirname(os.path.abspath(__file__))
sPath_library = os.path.dirname(os.path.dirname(sPath_current))
sys.path.append(sPath_library)

# Construct the relative path to the data folder
sFolder_data = os.path.join(sPath_current, '..', '..', 'data')
sFolder_data = os.path.realpath(sFolder_data)
# Print or use the data folder path
print(f"Data folder path: {sFolder_data}")


sFilename_mpas_mesh_in = os.path.join(sFolder_data, 'global','base_mesh.nc')

sFilename_png_out = os.path.join(sFolder_data, 'global', 'mpas_mesh.png')
# Test the culling fix with explicit parameters
visualize_mpas_mesh(sFilename_mpas_mesh_in,
                   sFilename_out = None,
                   iFlag_wireframe_only=True,
                   iFlag_cull_backfaces=True,
                   sCulling_mode='back',  # Use auto mode for best results
                   iFlag_verbose_in=True)  # Enable verbose output to see what's happening