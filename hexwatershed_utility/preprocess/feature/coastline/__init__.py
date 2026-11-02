"""
Coastline, land-ocean mask, and coastal flowline alignment utilities.
"""

from hexwatershed_utility.preprocess.feature.coastline.create_land_ocean_mask import (
    create_land_ocean_mask,
    create_land_ocean_mask_from_hydrobasins,
    create_land_ocean_mask_from_naturalearth,
    create_land_ocean_mask_from_shoreline,
    remove_internal_polygon,
)
from hexwatershed_utility.preprocess.feature.coastline.extend_land_mask_to_rivers import (
    extend_land_mask_to_rivers,
    fix_naturalearth_hydrosheds_incompatibility,
    geometries_bbox_overlap,
)

__all__ = [
    "create_land_ocean_mask",
    "create_land_ocean_mask_from_hydrobasins",
    "create_land_ocean_mask_from_naturalearth",
    "create_land_ocean_mask_from_shoreline",
    "extend_land_mask_to_rivers",
    "fix_naturalearth_hydrosheds_incompatibility",
    "geometries_bbox_overlap",
    "remove_internal_polygon",
]
