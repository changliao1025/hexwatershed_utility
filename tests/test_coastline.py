import os
import tempfile
import unittest

from osgeo import ogr, osr

from hexwatershed_utility.preprocess.feature.coastline import (
    create_land_ocean_mask,
    create_land_ocean_mask_from_basins,
    create_land_ocean_mask_from_naturalearth,
    create_land_ocean_mask_from_shoreline,
    extend_land_mask_to_rivers,
    fix_naturalearth_hydrosheds_incompatibility,
    geometries_bbox_overlap,
    remove_internal_polygon,
)


class TestCoastline(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_callables_exist(self):
        self.assertTrue(callable(create_land_ocean_mask_from_basins))
        self.assertTrue(callable(create_land_ocean_mask_from_shoreline))
        self.assertTrue(callable(create_land_ocean_mask_from_naturalearth))
        self.assertTrue(callable(create_land_ocean_mask))
        self.assertTrue(callable(extend_land_mask_to_rivers))
        self.assertTrue(callable(fix_naturalearth_hydrosheds_incompatibility))
        self.assertTrue(callable(geometries_bbox_overlap))
        self.assertTrue(callable(remove_internal_polygon))

    def test_geometries_bbox_overlap(self):
        bbox1 = (0.0, 10.0, 0.0, 10.0)
        bbox2 = (5.0, 15.0, 5.0, 15.0)
        bbox3 = (20.0, 30.0, 20.0, 30.0)
        self.assertTrue(geometries_bbox_overlap(bbox1, bbox2))
        self.assertFalse(geometries_bbox_overlap(bbox1, bbox3))

    def test_extend_land_mask_file_not_found(self):
        with self.assertRaises(FileNotFoundError):
            extend_land_mask_to_rivers([], "non_existent_land.geojson", "out.geojson")

    def test_create_land_ocean_mask_invalid_args(self):
        with self.assertRaises(ValueError):
            create_land_ocean_mask_from_basins("/tmp/out")

    def test_remove_internal_polygon(self):
        # Create a polygon with a hole
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)

        driver = ogr.GetDriverByName("GeoJSON")
        in_geojson = os.path.join(self.temp_dir.name, "with_hole.geojson")
        ds = driver.CreateDataSource(in_geojson)
        layer = ds.CreateLayer("land", srs, ogr.wkbPolygon)
        layer.CreateField(ogr.FieldDefn("name", ogr.OFTString))

        # 10x10 square with a 6x6 hole inside (area without hole = 100, with hole = 64)
        wkt_with_hole = "POLYGON ((0 0, 10 0, 10 10, 0 10, 0 0), (2 2, 8 2, 8 8, 2 8, 2 2))"
        geom_hole = ogr.CreateGeometryFromWkt(wkt_with_hole)
        feat = ogr.Feature(layer.GetLayerDefn())
        feat.SetGeometry(geom_hole)
        feat.SetField("name", "ContinentWithLake")
        layer.CreateFeature(feat)
        ds = None

        out_geojson = os.path.join(self.temp_dir.name, "without_hole.geojson")
        remove_internal_polygon(in_geojson, out_geojson)
        self.assertTrue(os.path.exists(out_geojson))

        ds_out = ogr.Open(out_geojson)
        layer_out = ds_out.GetLayer(0)
        self.assertEqual(layer_out.GetFeatureCount(), 1)
        feat_out = layer_out.GetFeature(0)
        geom_out = feat_out.GetGeometryRef()

        # Geometry count for a polygon is the number of rings (1 = exterior only)
        self.assertEqual(geom_out.GetGeometryCount(), 1)
        self.assertAlmostEqual(geom_out.GetArea(), 100.0, places=4)
        self.assertEqual(feat_out.GetField("name"), "ContinentWithLake")
        ds_out = None


if __name__ == "__main__":
    unittest.main()
