"""
River network feature utilities.
"""

from hexwatershed_utility.preprocess.feature.river.get_outlet_location import get_outlet_location
from hexwatershed_utility.preprocess.feature.river.simplify_river_network import (
    simplify_river_network,
)
from hexwatershed_utility.preprocess.feature.river.tag_river_outlet import tag_river_outlet

__all__ = [
    "get_outlet_location",
    "simplify_river_network",
    "tag_river_outlet",
]
