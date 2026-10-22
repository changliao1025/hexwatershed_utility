"""
Dam and reservoir feature utilities.
"""
__all__ = []

try:
    from hexwatershed_utility.preprocess.feature.dam.snap_dams_to_river_network import (
        snap_dams_to_river_networks,
    )
    __all__.append("snap_dams_to_river_networks")
except ImportError:
    pass
