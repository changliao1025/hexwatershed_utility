import os
from hexwatershed_utility.preprocess.mesh.visualize_mpas_mesh import visualize_mpas_mesh


sFilename_mpas_mesh_in = '/compyfs/liao313/04model/pyhexwatershed/global/pyflowline20251121001/jigsaw/out/base_mesh.nc'

sFilename_png_out = '/qfs/people/liao313/workspace/python/hexwatershed_utility/figures/mpas/mpas_mesh.png'
visualize_mpas_mesh(sFilename_mpas_mesh_in, sFilename_out = sFilename_png_out)