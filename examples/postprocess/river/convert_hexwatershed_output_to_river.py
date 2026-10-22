

import os

from hexwatershed_utility.postprocess.river.convert_hexwatershed_outputs_to_river_inputs import convert_hexwatershed_json_to_river_parameter_file

sFilename_json_in = '/compyfs/liao313/04model/pyhexwatershed/global/pyhexwatershed20260601005/hexwatershed/hexwatershed.json'

sFilename_river_parameter_in = '/compyfs/inputdata/rof/mosart/MOSART_Global_half_20210616.nc'

sWorkspace_out = '/compyfs/liao313/04model/pyhexwatershed/global/pyhexwatershed20260601005'

sFilename_river_parameter_out = os.path.join(sWorkspace_out, 'mosart_global_parameter_ocn30lnd10.nc')
sFilename_river_domain_out = os.path.join(sWorkspace_out, 'mosart_global_domain_ocn30lnd10.nc')




convert_hexwatershed_json_to_river_parameter_file(sFilename_json_in,
                                                   sFilename_river_parameter_in,
                                                   sFilename_river_parameter_out,
                                                   sFilename_river_domain_out)
