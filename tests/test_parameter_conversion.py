import os
import tempfile
import unittest

import netCDF4 as nc
import numpy as np

from hexwatershed_utility.postprocess.river.convert_hexwatershed_outputs_to_river_inputs import (
    convert_hexwatershed_json_to_river_parameter_file,
)
from tests.fixtures.synthetic_data import create_synthetic_mesh_json, create_synthetic_mosart_files


class TestParameterConversion(unittest.TestCase):
    def test_convert_hexwatershed_json_to_river_parameter_file(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            json_file = os.path.join(tmpdir, "hexwatershed.json")
            template_dom = os.path.join(tmpdir, "template_dom.nc")
            template_param = os.path.join(tmpdir, "template_param.nc")
            out_param = os.path.join(tmpdir, "river_param_out.nc")
            out_domain = os.path.join(tmpdir, "river_domain_out.nc")

            create_synthetic_mesh_json(json_file, n_cells=7)
            create_synthetic_mosart_files(template_dom, template_param, n_cells=7)

            p_out, d_out = convert_hexwatershed_json_to_river_parameter_file(
                sFilename_json_in=json_file,
                sFilename_river_parameter_in=template_param,
                sFilename_river_parameter_out=out_param,
                sFilename_river_domain_out=out_domain,
            )

            self.assertTrue(os.path.exists(p_out))
            self.assertTrue(os.path.exists(d_out))

            # Verify parameter file
            with nc.Dataset(p_out, "r") as ds:
                self.assertIn("gridcell", ds.dimensions)
                self.assertEqual(len(ds.dimensions["gridcell"]), 7)

                # Check critical variables
                expected_vars = [
                    "area",
                    "areaTotal",
                    "areaTotal2",
                    "dnID",
                    "CellID",
                    "ID",
                    "latixy",
                    "longxy",
                    "rdep",
                    "rlen",
                    "rwid",
                ]
                for var_name in expected_vars:
                    self.assertIn(var_name, ds.variables)

                # Cell 7 is outlet -> -9999 (masked in netCDF4)
                dn_id = ds.variables["dnID"][:]
                self.assertTrue(dn_id.data[6] == -9999 or np.ma.is_masked(dn_id[6]))

            # Verify domain file
            with nc.Dataset(d_out, "r") as ds_dom:
                self.assertTrue("ni" in ds_dom.dimensions or "gridcell" in ds_dom.dimensions)
                self.assertIn("xv", ds_dom.variables)
                self.assertIn("yv", ds_dom.variables)


if __name__ == "__main__":
    unittest.main()
