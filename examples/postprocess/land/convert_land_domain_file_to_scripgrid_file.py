import json
import os
from pathlib import Path

import netCDF4 as nc
import numpy as np
from pyesm.mesh.unstructured.esm_convert_unstructured_domain_file_to_scripgrid_file import (
    esm_convert_unstructured_domain_file_to_scripgrid_file,
)

# because elm domain only need to cell center and vertex info, it is best to use the mesh info file


sFilename_elm_domain_out = "/compyfs/liao313/04model/pyhexwatershed/global/pyhexwatershed20251208002/elm/elm_domain_from_hexwatershed.nc"

esm_convert_unstructured_domain_file_to_scripgrid_file(
    sFilename_elm_domain_out, sFilename_elm_domain_out.replace(".nc", ".scripgrid.nc")
)
