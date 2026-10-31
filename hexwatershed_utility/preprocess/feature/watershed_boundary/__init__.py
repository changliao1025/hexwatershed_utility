"""
Watershed boundary identification, extraction, and merging utilities.
"""

from hexwatershed_utility.preprocess.feature.watershed_boundary.find_boundary import (
    find_containing_boundary,
)
from hexwatershed_utility.preprocess.feature.watershed_boundary.merge_boundary import (
    merge_boundary_files,
    merge_boundary_polygons,
)

__all__ = [
    "find_containing_boundary",
    "merge_boundary_polygons",
    "merge_boundary_files",
]
