import subprocess
import os
import time
import numpy as np
import xarray as xr
from geometric_features import GeometricFeatures
from pyearth.system.python.retrieve_python_environment import retrieve_python_environment
from pyflowline.mesh.mpas.mpas_tools.mpasmsh import jigsaw_mesh_to_netcdf, inject_edge_tags, subtract_critical_passages, mask_reachable_ocean

import mpas_tools
from mpas_tools.mesh.conversion import convert, cull
from mpas_tools.logging import check_call
from mpas_tools.io import write_netcdf

HERE = os.path.abspath(os.path.dirname(__file__))
def debug_jigsaw_mpas_tools(sWorkspace_jigsaw_out, sFilename_triangles):
    netcdfFormat = "NETCDF4_CLASSIC" #NETCDF3_64BIT

    print("Forming base_mesh.nc")
    sFilename_base_mesh = os.path.join(sWorkspace_jigsaw_out, "out", "base_mesh_debug.nc")
    dummy = convert(xr.open_dataset( sFilename_triangles))
    write_netcdf(dummy, fileName=sFilename_base_mesh, format= netcdfFormat)

if __name__ == '__main__':
    #/compyfs/liao313/04model/pyflowline/arctic/pyflowline20250101001/jigsaw/tmp
    sWorkspace_jigsaw_out = '/compyfs/liao313/04model/pyflowline/arctic/pyflowline20250101001/jigsaw'
    sFilename_triangles = '/compyfs/liao313/04model/pyflowline/arctic/pyflowline20250101001/jigsaw/tmp/mesh_triangles.nc'
    debug_jigsaw_mpas_tools(sWorkspace_jigsaw_out, sFilename_triangles)
    pass
