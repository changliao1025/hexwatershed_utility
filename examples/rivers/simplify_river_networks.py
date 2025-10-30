import os
from hexwatershed_utility.preprocess.features.rivers.simplify_hydrorivers_networks import simplify_hydrorivers_networks
sFilename_flowline_hydroshed_in = '/compyfs/liao313/00raw/hydrology/hydrosheds/hydroriver/HydroRIVERS_v10_na_shp/HydroRIVERS_v10_na.shp'

#global
#sFilename_flowline_hydroshed_in = '/compyfs/liao313/00raw/hydrology/hydrosheds/hydroriver/HydroRIVERS_v10_shp/HydroRIVERS_v10_shp/HydroRIVERS_v10.shp'

sWorkspace_out = '/qfs/people/liao313/data/hexwatershed/conus/vector/hydrology'

sWorkspace_out = '/compyfs/liao313/00raw/hydrology/conus/vector/rivers'

dResolution_land = 10
dDistance_tolerance_in = dResolution_land * 1.0E3


dDrainage_area_threshold_in = dResolution_land * dResolution_land *10 * 1.0E6 #km2
sDistance_tolerance = "{:.2E}".format(dDistance_tolerance_in)
sDrainage_area_threshold = "{:.2E}".format(dDrainage_area_threshold_in)

sFoldname = sDistance_tolerance + '_' + sDrainage_area_threshold

sFilename_flowline_hydroshed_out = 'HydroRIVERS_v10_' + sDistance_tolerance + '_' + sDrainage_area_threshold + '_simplified.geojson'

sWorkspace_out = os.path.join(sWorkspace_out, sFoldname)
if not os.path.exists(sWorkspace_out):
    os.makedirs(sWorkspace_out)

sFilename_flowline_hydroshed_out =  os.path.join(sWorkspace_out, sFilename_flowline_hydroshed_out)

simplify_hydrorivers_networks(sFilename_flowline_hydroshed_in, sFilename_flowline_hydroshed_out, dDistance_tolerance_in, dDrainage_area_threshold_in)
