"""
Mesh pre-processing and format conversion utilities.
"""

__all__ = []

try:
    from hexwatershed_utility.preprocess.mesh.convert_mpas_mesh_to_vector import (
        convert_mpas_mesh_to_vector,
    )

    __all__.append("convert_mpas_mesh_to_vector")
except ImportError:
    pass

try:
    from hexwatershed_utility.preprocess.mesh.inject_dem_into_mpas_mesh import (
        inject_dem_into_mpas_mesh_netcdf,
    )

    __all__.append("inject_dem_into_mpas_mesh_netcdf")
except ImportError:
    pass
