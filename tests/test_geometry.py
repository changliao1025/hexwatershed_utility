import unittest

import numpy as np

from hexwatershed_utility.postprocess.river.get_geometry import get_geometry


class TestGeometry(unittest.TestCase):
    def test_get_geometry_empirical_mode(self):
        # 3 cells
        aLongitude = np.array([-100.0, -99.5, -99.0])
        aLatitude = np.array([40.0, 40.5, 41.0])
        aCellID = np.array([1, 2, 3])
        aCellID_downslope = np.array([2, 3, -1])
        aArea = np.array([1.0e6, 2.0e6, 3.0e6])  # m²

        aWidth, aDepth, aFlood_2yr = get_geometry(
            aLongitude_in=aLongitude,
            aLatitude_in=aLatitude,
            aCellID=aCellID,
            aCellID_downslope=aCellID_downslope,
            aArea=aArea,
            iFlag_read_runoff=0,
        )

        self.assertEqual(len(aWidth), 3)
        self.assertEqual(len(aDepth), 3)
        self.assertEqual(len(aFlood_2yr), 3)

        # All values should be positive
        self.assertTrue(np.all(aWidth > 0))
        self.assertTrue(np.all(aDepth > 0))
        self.assertTrue(np.all(aFlood_2yr > 0))

        # Larger area should yield larger flood discharge, width, depth
        self.assertTrue(aFlood_2yr[0] < aFlood_2yr[1] < aFlood_2yr[2])
        self.assertTrue(aWidth[0] < aWidth[1] < aWidth[2])
        self.assertTrue(aDepth[0] < aDepth[1] < aDepth[2])


if __name__ == "__main__":
    unittest.main()
