import os

from hexwatershed_utility.preprocess.feature.river.simplify_hydrorivers_network import (
    simplify_hydrorivers_network,
)

sFilename_flowline_hydroshed_in = "/compyfs/liao313/00raw/hydrology/hydrosheds/hydroriver/HydroRIVERS_v10_na_shp/HydroRIVERS_v10_na.shp"

# global
sFilename_flowline_hydroshed_in = "/compyfs/liao313/00raw/hydrology/hydrosheds/hydroriver/HydroRIVERS_v10_shp/HydroRIVERS_v10_shp/HydroRIVERS_v10.shp"

sWorkspace_out = "/qfs/people/liao313/data/hexwatershed/conus/vector/hydrology"

sWorkspace_out = "/compyfs/liao313/00raw/hydrology/global/vector/rivers"

dResolution_land = 10
dDistance_tolerance_in = dResolution_land * 1.0e3  # in meter


dDrainage_area_threshold_in = dResolution_land * dResolution_land * 100 * 1.0e6  # m2
sDistance_tolerance = f"{dDistance_tolerance_in:.2E}"
sDrainage_area_threshold = f"{dDrainage_area_threshold_in:.2E}"

sFoldname = sDistance_tolerance + "_" + sDrainage_area_threshold

sFilename_flowline_hydroshed_out = (
    "HydroRIVERS_v10_"
    + sDistance_tolerance
    + "_"
    + sDrainage_area_threshold
    + "_simplified.geojson"
)

sWorkspace_out = os.path.join(sWorkspace_out, sFoldname)
if not os.path.exists(sWorkspace_out):
    os.makedirs(sWorkspace_out)

sFilename_flowline_hydroshed_out = os.path.join(
    sWorkspace_out, sFilename_flowline_hydroshed_out
)

simplify_hydrorivers_networks(
    sFilename_flowline_hydroshed_in,
    sFilename_flowline_hydroshed_out,
    dDistance_tolerance_in,
    dDrainage_area_threshold_in,
    iFlag_pyflowline_configuration_in=0,
    nOutlet_largest=100,
)
