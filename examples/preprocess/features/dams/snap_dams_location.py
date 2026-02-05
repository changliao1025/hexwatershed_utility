from hexwatershed_utility.preprocess.features.dams.snap_dams_to_river_networks import snap_dams_to_river_networks

sWorkspace_data = '/compyfs/liao313/00raw/dam/GOODD_data/Data'



sFilename_river= '/qfs/people/liao313/data/hexwatershed/global/vector/river_networks_wo_lakes.geojson'

sFilename_dam="/compyfs/zhou014/datasets/Reservoir_database/GRanD_Version_1_3/GRanD_dams_v1_3.shp"
sFilename_river="/compyfs/liao313/04model/pyhexwatershed/global/river_network/1.00E+04_1.00E+10/HydroRIVERS_v10_simplified_1.00E+04_1.00E+10.geojson"
sFilename_out="/compyfs/liao313/00raw/dam/GRanD_Version_1_3/GRanD_dams_v1_3.geojson"

sVar_drainage='drainage' #Catch_skm
sVar_drainage = 'Catch_skm'

# The function will create three output files:
# 1. GRanD_dams_v1_3_snapped.geojson - Dams successfully snapped to river networks (distance < 3km)
# 2. GRanD_dams_v1_3_no_river.geojson - Dams not close to any river
# 3. GRanD_dams_v1_3_too_far.geojson - Dams too far from rivers (distance >= 3km)

snap_dams_to_river_networks(sFilename_dam, sFilename_river, sFilename_out,
                             sVariable_drainage=sVar_drainage)