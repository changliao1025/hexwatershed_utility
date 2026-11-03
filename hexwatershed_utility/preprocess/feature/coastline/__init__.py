"""
Coastline, land-ocean mask, and coastal flowline alignment utilities.
"""

from hexwatershed_utility.preprocess.feature.coastline.create_land_ocean_mask import (
    create_land_ocean_mask_from_hydrobasins,
    create_land_ocean_mask_from_naturalearth,
    simplify_vector_geometry,
)
from hexwatershed_utility.preprocess.feature.coastline.extend_land_mask_to_rivers import (
    extend_land_mask_to_rivers,
    geometries_bbox_overlap,
)

__all__ = [
    "create_land_ocean_mask_from_hydrobasins",
    "create_land_ocean_mask_from_naturalearth",
    "extend_land_mask_to_rivers",
    "geometries_bbox_overlap",
    "simplify_vector_geometry",
]
