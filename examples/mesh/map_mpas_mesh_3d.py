import os
from hexwatershed_utility.preprocess.mesh.map_mpas_mesh_on_sphere import map_mpas_mesh_on_sphere


sFilename_mpas_mesh_in = os.path.join('/compyfs/liao313/04model/pyhexwatershed/global/pyflowline20260203001/jigsaw/out','base_mesh.nc')
sFilename_png_out = os.path.join('/qfs/people/liao313/workspace/python/hexwatershed_utility/figures', 'mpas_mesh_ocn6_18coast6lnd12riv6_3d.jpg')

# Test the culling fix with explicit parameters
map_mpas_mesh_on_sphere(sFilename_mpas_mesh_in,
                   sFilename_out = sFilename_png_out,
                   base_layer = 'natural_earth_1',
                    dLongitude_focus_in= -135,
                dLatitude_focus_in= 45,
                   iFlag_wireframe_only = True,
                  window_size_in=(8000, 8000),  # Set the window size for better visibility
                   iFlag_verbose_in=True)  # Enable verbose output to see what's happening