from hexwatershed_utility.preprocess.feature.watershed_boundary import merge_boundary_files

# use this example to merge multiple subbasin as a single watershed boundary

aFilename = list()

# add them one by one

aFilename.append("/Users/changliao/data/hexwatershed/yellowriver/vector/basin_boundary/01.geojson")
aFilename.append("/Users/changliao/data/hexwatershed/yellowriver/vector/basin_boundary/02.geojson")
aFilename.append("/Users/changliao/data/hexwatershed/yellowriver/vector/basin_boundary/03.geojson")
aFilename.append("/Users/changliao/data/hexwatershed/yellowriver/vector/basin_boundary/04.geojson")
aFilename.append("/Users/changliao/data/hexwatershed/yellowriver/vector/basin_boundary/05.geojson")


merge_boundary_files(
    aFilename, "/Users/changliao/data/hexwatershed/yellowriver/vector/basin_boundary/merged.geojson"
)
