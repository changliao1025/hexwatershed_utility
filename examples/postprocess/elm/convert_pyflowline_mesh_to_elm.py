import os
from pathlib import Path
import json
import numpy as np
import netCDF4 as nc
from pye3sm.mesh.unstructured.e3sm_create_unstructured_domain_file_full import e3sm_create_unstructured_domain_file_full
from pye3sm.mesh.unstructured.e3sm_convert_unstructured_domain_file_to_scripgrid_file import e3sm_convert_unstructured_domain_file_to_scripgrid_file

#because elm domain only need to cell center and vertex info, it is best to use the mesh info file

from hexwatershed_utility.postprocess.elm.convert_hexwatershed_outputs_elm_inputs import convert_hexwatershed_json_to_elm_domain_file
sFilename_json_in = '/compyfs/liao313/04model/pyhexwatershed/global/pyhexwatershed20251208002/hexwatershed/hexwatershed.json'
sFilename_elm_domain_out = '/compyfs/liao313/04model/pyhexwatershed/global/pyhexwatershed20251208002/elm/elm_domain_from_hexwatershed.nc'

convert_hexwatershed_json_to_elm_domain_file(sFilename_json_in,
                                               sFilename_elm_domain_out)