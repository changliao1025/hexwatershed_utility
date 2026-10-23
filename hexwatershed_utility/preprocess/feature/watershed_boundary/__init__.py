"""
Watershed boundary extraction and repair utilities.
"""

__all__ = []

try:
    from hexwatershed_utility.preprocess.feature.watershed_boundary.reorganize_boundary import (
        reorganize_boundary,
    )

    __all__.append("reorganize_boundary")
except ImportError:
    pass

try:
    from hexwatershed_utility.preprocess.feature.watershed_boundary.shapefile_repair_utility import (
        diagnose_and_fix_shapefile,
        open_shapefile_with_fallback,
    )

    __all__.extend(["open_shapefile_with_fallback", "diagnose_and_fix_shapefile"])
except ImportError:
    pass

try:
    from hexwatershed_utility.preprocess.feature.watershed_boundary.find_minimal_hydrobasins_watershed_boundary import (
        find_minimal_hydrobasins_watershed_boundary,
    )

    __all__.append("find_minimal_hydrobasins_watershed_boundary")
except ImportError:
    pass
