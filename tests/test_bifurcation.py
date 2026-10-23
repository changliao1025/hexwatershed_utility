import os
import tempfile
import unittest

try:
    import pandas as pd  # noqa: F401

    from hexwatershed_utility.postprocess.river.create_bifurc_netcdf import (
        create_bifurc_netcdf,
        verify_bifurc_file,
    )

    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False

import netCDF4 as nc

from tests.fixtures.synthetic_data import (
    create_synthetic_bifurcation_csv,
    create_synthetic_mosart_files,
)


@unittest.skipUnless(HAS_PANDAS, "pandas is required for bifurcation NetCDF creation")
class TestBifurcation(unittest.TestCase):
    def test_create_and_verify_bifurc_netcdf(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            dom_nc = os.path.join(tmpdir, "domain.nc")
            param_nc = os.path.join(tmpdir, "param.nc")
            bifurc_nc = os.path.join(tmpdir, "bifurc.nc")
            csv_path = os.path.join(tmpdir, "bifurc.csv")

            create_synthetic_mosart_files(dom_nc, param_nc, n_cells=7)
            create_synthetic_bifurcation_csv(csv_path)

            create_bifurc_netcdf(
                input_file=param_nc,
                output_file=bifurc_nc,
                csv_file=csv_path,
                max_downstream=2,
            )

            self.assertTrue(os.path.exists(bifurc_nc))

            # Inspect NetCDF variables and dimensions
            with nc.Dataset(bifurc_nc, "r") as ds:
                self.assertIn("downstream", ds.dimensions)
                self.assertEqual(len(ds.dimensions["downstream"]), 2)
                self.assertIn("dnID_original", ds.variables)
                self.assertIn("bifurc_ratio", ds.variables)
                self.assertIn("ibt_demand", ds.variables)

                ratio = ds.variables["bifurc_ratio"][:]
                self.assertEqual(ratio.shape, (7, 2))

            # Verify with verification utility
            verify_bifurc_file(bifurc_nc, num_bifurc_points=1, is_structured=False)


if __name__ == "__main__":
    unittest.main()
