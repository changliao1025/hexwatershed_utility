import json
import os
from pathlib import Path

import netCDF4 as nc
import numpy as np
from pyesm.mesh.unstructured.esm_convert_unstructured_domain_file_to_scripgrid_file import (
    esm_convert_unstructured_domain_file_to_scripgrid_file,
)
from pyesm.mesh.unstructured.esm_create_unstructured_domain_file_full import (
    esm_create_unstructured_domain_file_full,
)

from hexwatershed_utility.postprocess.land.convert_hexwatershed_outputs_land_inputs import (
    convert_hexwatershed_json_to_land_domain_file,
)

# because elm domain only need to cell center and vertex info, it is best to use the mesh info file


sFilename_json_in = "/compyfs/liao313/04model/pyhexwatershed/global/pyhexwatershed20251208002/hexwatershed/hexwatershed.json"
sFilename_land_domain_out = "/compyfs/liao313/04model/pyhexwatershed/global/pyhexwatershed20251208002/elm/elm_domain_from_hexwatershed.nc"

convert_hexwatershed_json_to_land_domain_file(
    sFilename_json_in, sFilename_land_domain_out
)
