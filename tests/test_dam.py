import unittest

from hexwatershed_utility.preprocess.feature.dam import snap_dams_to_river_network


class TestDam(unittest.TestCase):
    def test_file_not_found(self):
        with self.assertRaises(FileNotFoundError):
            snap_dams_to_river_network(
                "non_existent_dam.shp", "non_existent_river.geojson", "out.geojson"
            )


if __name__ == "__main__":
    unittest.main()
