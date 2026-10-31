"""
Example: Find containing watershed boundary for a river network.
"""

from hexwatershed_utility.preprocess.feature.watershed_boundary import (
    find_containing_boundary,
    merge_boundary_polygons,
)

# 1. Input paths
sFilename_river_network = "/path/to/river_network.geojson"
sFolder_boundary = "/path/to/hydrobasins"  # or direct list of boundary files
sFilename_out = "/path/to/output_watershed_boundary.geojson"

# 2. Find containing boundary (single union boundary or list of subbasin boundaries)
wkt_result = find_containing_boundary(
    sFilename_river_network,
    sFolder_boundary_in=sFolder_boundary,
    iFlag_nested_in=False,  # Set True to extract per-subbasin boundaries
)

# 3. Save to GeoJSON (merge_boundary_polygons handles deduplication & export in 1 line)
if wkt_result is not None:
    wkt_list = [wkt_result] if isinstance(wkt_result, str) else wkt_result
    merge_boundary_polygons(wkt_list, sFilename_out=sFilename_out)
    print(f"Watershed boundary saved to: {sFilename_out}")
else:
    print("No containing boundary found.")
