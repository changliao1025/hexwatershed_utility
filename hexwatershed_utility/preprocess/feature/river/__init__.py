"""
River network feature utilities.
"""

__all__ = []

try:
    from hexwatershed_utility.preprocess.feature.river.get_outlet_location import (
        get_outlet_location,
    )

    __all__.append("get_outlet_location")
except ImportError:
    pass

try:
    from hexwatershed_utility.preprocess.feature.river.simplify_hydrorivers_network import (
        simplify_hydrorivers_network,
    )

    __all__.append("simplify_hydrorivers_network")
except ImportError:
    pass

try:
    from hexwatershed_utility.preprocess.feature.river.tag_river_outlet import (
        tag_river_outlet,
    )

    __all__.append("tag_river_outlet")
except ImportError:
    pass
