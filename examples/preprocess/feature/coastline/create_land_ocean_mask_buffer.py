"""
Example demonstrating how to generate a land-ocean mask and coastal buffer
from watershed basin boundaries (e.g., HydroBASINS or regional basins).
"""

import os

from hexwatershed_utility.preprocess.feature.coastline import create_land_ocean_mask_from_basins


def main():
    sWorkspace_output = os.path.expanduser("~/tmp/coastline_basin_mask_example")
    os.makedirs(sWorkspace_output, exist_ok=True)

    # Example: Generate land mask and coastal buffer zone from basin polygons
    # Can accept an explicit list of files or a folder with hierarchy level
    sFolder_basins = os.path.expanduser("~/data/hydrobasins")

    if os.path.exists(sFolder_basins):
        sFilename_tif, sFilename_buffer_parquet = create_land_ocean_mask_from_basins(
            sWorkspace_output=sWorkspace_output,
            sFolder_boundary_in=sFolder_basins,
            sLevel="lev03",
            dResolution_x=30 / 3600.0,
            dResolution_y=30 / 3600.0,
            dThreshold_area_island=1000.0,
            dResolution_coastline_buffer=10.0,  # 10 km buffer using pyearthbuffer
        )
        print("Basin mask and buffer completed:")
        print(f"  Raster GeoTIFF   : {sFilename_tif}")
        print(f"  Buffered Parquet : {sFilename_buffer_parquet}")
    else:
        print(
            f"Sample basin directory {sFolder_basins} not found. Please provide valid basin boundary files."
        )


if __name__ == "__main__":
    main()
