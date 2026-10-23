import unittest

from osgeo import ogr

from hexwatershed_utility.preprocess.feature.watershed_boundary.reorganize_boundary import (
    reorganize_boundary,
)


class TestBoundary(unittest.TestCase):
    def test_reorganize_boundary_empty(self):
        wkt = reorganize_boundary([])
        self.assertIn("MULTIPOLYGON", wkt.upper())

    def test_reorganize_boundary_single_polygon(self):
        poly_wkt = "POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0))"
        wkt = reorganize_boundary([poly_wkt])
        self.assertIn("MULTIPOLYGON", wkt.upper())
        geom = ogr.CreateGeometryFromWkt(wkt)
        self.assertIsNotNone(geom)
        self.assertAlmostEqual(geom.GetArea(), 1.0, places=5)

    def test_reorganize_boundary_duplicate_polygons(self):
        # Two identical polygons - one should be deduplicated
        poly1 = "POLYGON ((0 0, 2 0, 2 2, 0 2, 0 0))"
        poly2 = "POLYGON ((0 0, 2 0, 2 2, 0 2, 0 0))"
        wkt = reorganize_boundary([poly1, poly2])
        geom = ogr.CreateGeometryFromWkt(wkt)
        self.assertIsNotNone(geom)
        self.assertAlmostEqual(geom.GetArea(), 4.0, places=5)
        self.assertEqual(geom.GetGeometryCount(), 1)


if __name__ == "__main__":
    unittest.main()
