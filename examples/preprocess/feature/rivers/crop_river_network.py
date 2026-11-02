from hexwatershed_utility.preprocess.feature.river.crop_river_network_by_boundary import (
    crop_river_network_by_boundary,
)

sFilename_river = "/Users/changliao/scratch/rawdata/hydrology/hydrosheds/hydroriver/HydroRIVERS_v10_as_shp/HydroRIVERS_v10_as_shp/HydroRIVERS_v10_as.shp"
sFilename_boundary = (
    "/Users/changliao/data/hexwatershed/yellowriver/vector/basin_boundary/merged.geojson"
)
sFilename_out = (
    "/Users/changliao/data/hexwatershed/yellowriver/vector/river_network/cropped_river_network.shp"
)

crop_river_network_by_boundary(
    sFilename_river=sFilename_river,
    sFilename_boundary=sFilename_boundary,
    sFilename_out=sFilename_out,
)
