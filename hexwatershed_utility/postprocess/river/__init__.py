"""
River transport model post-processing utilities.
"""
from hexwatershed_utility.postprocess.river.get_geometry import get_geometry
from hexwatershed_utility.postprocess.river.find_contributing_cells import (
    find_contributing_cells,
    build_contributing_cells_all,
)

__all__ = [
    "get_geometry",
    "find_contributing_cells",
    "build_contributing_cells_all",
]

try:
    from hexwatershed_utility.postprocess.river.convert_hexwatershed_outputs_to_river_inputs import (
        convert_hexwatershed_json_to_river_parameter_file,
    )
    __all__.append("convert_hexwatershed_json_to_river_parameter_file")
except ImportError:
    pass

try:
    from hexwatershed_utility.postprocess.river.create_bifurc_netcdf import (
        create_bifurc_netcdf,
    )
    __all__.append("create_bifurc_netcdf")
except ImportError:
    pass
