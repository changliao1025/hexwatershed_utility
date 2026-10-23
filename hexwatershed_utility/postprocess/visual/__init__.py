"""
Flow direction and post-processing visualization utilities.
"""

__all__ = []

try:
    from hexwatershed_utility.postprocess.visual.map_flow_direction import (
        map_flow_direction,
    )

    __all__.append("map_flow_direction")
except ImportError:
    pass

try:
    from hexwatershed_utility.postprocess.visual.map_flow_direction_on_sphere import (
        map_flow_direction_on_sphere,
    )

    __all__.append("map_flow_direction_on_sphere")
except ImportError:
    pass
