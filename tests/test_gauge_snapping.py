import os
import tempfile
import unittest

from hexwatershed_utility.postprocess.river.find_gage_mesh_cell_id import find_gage_mesh_cell_id
from tests.fixtures.synthetic_data import create_synthetic_mosart_files


class TestGaugeSnapping(unittest.TestCase):
    def test_find_gage_mesh_cell_id(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            dom_nc = os.path.join(tmpdir, "domain.nc")
            param_nc = os.path.join(tmpdir, "param.nc")
            create_synthetic_mosart_files(dom_nc, param_nc, n_cells=7)

            # Gauge near Cell 3 (-99.8, 40.0) with area 50 km2
            site_names = ["GAGE_03"]
            lon_gage = [-99.8]
            lat_gage = [40.0]
            drainage_area = [50.0]  # km2

            gage_idx, cell_ids, model_areas = find_gage_mesh_cell_id(
                aSitename_in=site_names,
                aLongitude_gage_in=lon_gage,
                aLatitude_gage_in=lat_gage,
                aDrainage_area_in=drainage_area,
                sFilename_domain_in=dom_nc,
                sFilename_parameter_in=param_nc,
                dThreshold_drainage_in=5.0e6,
                dThreshold_difference_in=0.15,
                iFlag_data_km_in=0,
                dBuffer_in=0.2,
            )

            self.assertEqual(len(gage_idx), 1)
            self.assertEqual(gage_idx[0], 0)
            self.assertEqual(cell_ids[0], 3)
            self.assertAlmostEqual(model_areas[0], 50.0, places=1)


if __name__ == "__main__":
    unittest.main()
