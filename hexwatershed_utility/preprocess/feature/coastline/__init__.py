"""
Coastline and land-ocean mask utilities.
"""
__all__ = []

try:
    from hexwatershed_utility.preprocess.feature.coastline.create_land_ocean_mask_from_hydrobasin import (
        create_land_ocean_mask_from_hydrobasin,
    )
    __all__.append("create_land_ocean_mask_from_hydrobasin")
except ImportError:
    pass

try:
    from hexwatershed_utility.preprocess.feature.coastline.create_land_ocean_mask_from_naturalearth import (
        create_land_ocean_mask_from_naturalearth,
    )
    __all__.append("create_land_ocean_mask_from_naturalearth")
except ImportError:
    pass

try:
    from hexwatershed_utility.preprocess.feature.coastline.fix_naturalearth_hydrosheds_incompatibility import (
        fix_naturalearth_hydrosheds_incompatibility,
    )
    __all__.append("fix_naturalearth_hydrosheds_incompatibility")
except ImportError:
    pass
