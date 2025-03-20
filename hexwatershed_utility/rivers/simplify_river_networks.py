from pyflowline.algorithms.simplification.simplify_hydroshed import simplify_hydroshed
sFilename_flowline_hydroshed_in = '/compyfs/liao313/00raw/hydrology/hydroshed/hydroriver/HydroRIVERS_v10_na_shp/HydroRIVERS_v10_na.shp'

#global
sFilename_flowline_hydroshed_in = '/compyfs/liao313/00raw/hydrology/hydroshed/hydroriver/HydroRIVERS_v10_shp/HydroRIVERS_v10_shp/HydroRIVERS_v10.shp'

sFilename_flowline_hydroshed_out = '/compyfs/liao313/04model/pyflowline/global/flowline_hydroshed_simplified_4.0E4.geojson'

dDistance_tolerance_in = 4.0 * 1.0E3
dDrainage_area_threshold_in = 4.0E4  * 1.0E6
simplify_hydroshed(sFilename_flowline_hydroshed_in, sFilename_flowline_hydroshed_out, dDistance_tolerance_in, dDrainage_area_threshold_in)
