import unittest

from osgeo import ogr

from hexwatershed_utility.preprocess.feature.watershed_boundary import (
    find_containing_boundary,
    merge_boundary_files,
    merge_boundary_polygons,
)
from hexwatershed_utility.preprocess.feature.watershed_boundary.find_boundary import (
    _find_boundary_for_geometry,
    _find_containing_feature_in_file,
)
from hexwatershed_utility.preprocess.feature.watershed_boundary.merge_boundary import (
    _extract_unique_polygons,
)


class TestBoundary(unittest.TestCase):
    def test_merge_boundary_polygons_empty(self):
        wkt = merge_boundary_polygons([])
        self.assertIn("MULTIPOLYGON", wkt.upper())

    def test_merge_boundary_polygons_single_polygon(self):
        poly_wkt = "POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0))"
        wkt = merge_boundary_polygons([poly_wkt])
        self.assertIn("MULTIPOLYGON", wkt.upper())
        geom = ogr.CreateGeometryFromWkt(wkt)
        self.assertIsNotNone(geom)
        self.assertAlmostEqual(geom.GetArea(), 1.0, places=5)

    def test_merge_boundary_polygons_duplicate_polygons(self):
        # Two identical polygons - one should be deduplicated
        poly1 = "POLYGON ((0 0, 2 0, 2 2, 0 2, 0 0))"
        poly2 = "POLYGON ((0 0, 2 0, 2 2, 0 2, 0 0))"
        wkt = merge_boundary_polygons([poly1, poly2])
        geom = ogr.CreateGeometryFromWkt(wkt)
        self.assertIsNotNone(geom)
        self.assertAlmostEqual(geom.GetArea(), 4.0, places=5)
        self.assertEqual(geom.GetGeometryCount(), 1)

    def test_merge_boundary_polygons_union(self):
        # Two adjacent squares touching along x=1 should union into a single rectangle
        poly1 = "POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0))"
        poly2 = "POLYGON ((1 0, 2 0, 2 1, 1 1, 1 0))"
        wkt = merge_boundary_polygons([poly1, poly2], iFlag_union=True)
        geom = ogr.CreateGeometryFromWkt(wkt)
        self.assertIsNotNone(geom)
        self.assertAlmostEqual(geom.GetArea(), 2.0, places=5)
        self.assertEqual(geom.GetGeometryCount(), 1)

    def test_extract_unique_polygons_empty(self):
        res = _extract_unique_polygons([])
        self.assertEqual(res, [])

    def test_find_containing_boundary_callable(self):
        self.assertTrue(callable(find_containing_boundary))
        self.assertTrue(callable(merge_boundary_polygons))
        self.assertTrue(callable(merge_boundary_files))


if __name__ == "__main__":
    unittest.main()
