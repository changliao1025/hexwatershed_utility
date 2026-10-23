import unittest

import numpy as np

from hexwatershed_utility.postprocess.river.calculate_metrics import (
    calculate_kge,
    calculate_nse,
    calculate_pbias,
    calculate_rmse,
)


class TestMetrics(unittest.TestCase):
    def setUp(self):
        self.obs = np.array([10.0, 20.0, 30.0, 40.0, 50.0])
        self.sim = np.array([10.0, 20.0, 30.0, 40.0, 50.0])

    def test_nse_perfect_match(self):
        nse = calculate_nse(self.obs, self.sim)
        self.assertAlmostEqual(nse, 1.0, places=5)

    def test_nse_constant_observation(self):
        obs = np.array([5.0, 5.0, 5.0])
        sim = np.array([5.0, 6.0, 5.0])
        self.assertTrue(np.isnan(calculate_nse(obs, sim)))

    def test_nse_with_nans(self):
        obs = np.array([10.0, np.nan, 30.0, 40.0, 50.0])
        sim = np.array([10.0, 20.0, 30.0, np.nan, 50.0])
        nse = calculate_nse(obs, sim)
        self.assertAlmostEqual(nse, 1.0, places=5)

    def test_kge_perfect_match(self):
        kge = calculate_kge(self.obs, self.sim)
        self.assertAlmostEqual(kge, 1.0, places=5)

    def test_kge_insufficient_data(self):
        self.assertTrue(np.isnan(calculate_kge([5.0], [5.0])))
        self.assertTrue(np.isnan(calculate_kge([], [])))

    def test_kge_zero_mean_or_std(self):
        obs = [0.0, 0.0, 0.0]
        sim = [1.0, 2.0, 3.0]
        self.assertTrue(np.isnan(calculate_kge(obs, sim)))

    def test_rmse_perfect_match(self):
        rmse = calculate_rmse(self.obs, self.sim)
        self.assertAlmostEqual(rmse, 0.0, places=5)

    def test_rmse_known_difference(self):
        sim = self.obs + 2.0
        rmse = calculate_rmse(self.obs, sim)
        self.assertAlmostEqual(rmse, 2.0, places=5)

    def test_pbias_perfect_match(self):
        pbias = calculate_pbias(self.obs, self.sim)
        self.assertAlmostEqual(pbias, 0.0, places=5)

    def test_pbias_overestimation(self):
        # sum(obs) = 150, sum(sim-obs) = 15 -> 10%
        sim = self.obs + 3.0
        pbias = calculate_pbias(self.obs, sim)
        self.assertAlmostEqual(pbias, 10.0, places=5)

    def test_empty_arrays(self):
        self.assertTrue(np.isnan(calculate_nse([], [])))
        self.assertTrue(np.isnan(calculate_rmse([], [])))
        self.assertTrue(np.isnan(calculate_pbias([], [])))


if __name__ == "__main__":
    unittest.main()
