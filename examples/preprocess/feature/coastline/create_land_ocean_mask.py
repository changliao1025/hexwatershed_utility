"""
Example demonstrating how to generate a global land-ocean mask from
shoreline vector data (e.g. Natural Earth) or basin boundaries (e.g. HydroBASINS).
"""

import os

from hexwatershed_utility.preprocess.feature.coastline import create_land_ocean_mask_from_shoreline


def main():
    sWorkspace_output = os.path.expanduser("/Users/changliao/data/hexwatershed/global/coastline")
    os.makedirs(sWorkspace_output, exist_ok=True)

    # Example 1: Generate land-ocean mask from shoreline vector data
    # (e.g., Natural Earth 10m land). By default, internal rings (inland lakes)
    # are preserved (iFlag_remove_internal_rings=False).
    sFilename_mask_tif, sFilename_mask_vector = create_land_ocean_mask_from_shoreline(
        sWorkspace_output=sWorkspace_output,
        sFilename_land_vector_in=None,  # Auto-downloads Natural Earth 10m land if None
        dResolution_x=0.1,  # 0.1 degree resolution for demo
        dResolution_y=0.1,
        dThreshold_area_island=1000.0,  # Filter out islands < 1000 km²
        iRaster_buffer_pixel=1,  # 1-pixel raster coastline buffer
        iFlag_remove_internal_rings=False,  # Retain interior water bodies by default
    )
    print("Shoreline mask completed:")
    print(f"  Raster GeoTIFF : {sFilename_mask_tif}")
    print(f"  Vector GeoJSON : {sFilename_mask_vector}")


if __name__ == "__main__":
    main()
