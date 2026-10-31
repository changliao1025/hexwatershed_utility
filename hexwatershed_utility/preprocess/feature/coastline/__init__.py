"""
Coastline, land-ocean mask, and coastal flowline alignment utilities.
"""

from hexwatershed_utility.preprocess.feature.coastline.create_land_ocean_mask import (
    create_land_ocean_mask_from_basins,
    create_land_ocean_mask_from_shoreline,
    remove_internal_polygon,
)
from hexwatershed_utility.preprocess.feature.coastline.extend_land_mask_to_rivers import (
    extend_land_mask_to_rivers,
)

__all__ = [
    "create_land_ocean_mask_from_basins",
    "create_land_ocean_mask_from_shoreline",
    "remove_internal_polygon",
    "extend_land_mask_to_rivers",
]
