import os
import tempfile
import unittest

from osgeo import ogr, osr

from hexwatershed_utility.preprocess.feature.river import (
    crop_river_network_by_boundary,
    get_outlet_location,
    simplify_river_network,
    tag_river_outlet,
)


class TestRiver(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_callables_exist(self):
        self.assertTrue(callable(crop_river_network_by_boundary))
        self.assertTrue(callable(get_outlet_location))
        self.assertTrue(callable(simplify_river_network))
        self.assertTrue(callable(tag_river_outlet))

    def test_crop_river_network_by_boundary_file_not_found(self):
        with self.assertRaises(FileNotFoundError):
            crop_river_network_by_boundary(
                "non_existent_river.geojson", "boundary.geojson", "out.geojson"
            )

    def test_get_outlet_location_file_not_found(self):
        with self.assertRaises(FileNotFoundError):
            get_outlet_location("non_existent_flowline.geojson")

    def test_get_outlet_location_basic(self):
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)

        driver = ogr.GetDriverByName("GeoJSON")
        shp_path = os.path.join(self.temp_dir.name, "test_river.geojson")
        ds = driver.CreateDataSource(shp_path)
        layer = ds.CreateLayer("river", srs, ogr.wkbLineString)
        layer.CreateField(ogr.FieldDefn("downstream_id", ogr.OFTInteger))

        # Flowline 1 (tributary)
        line1 = ogr.CreateGeometryFromWkt("LINESTRING (-80 40, -80.5 40.5)")
        feat1 = ogr.Feature(layer.GetLayerDefn())
        feat1.SetGeometry(line1)
        feat1.SetField("downstream_id", 2)
        layer.CreateFeature(feat1)

        # Flowline 2 (outlet reach, downstream_id = 0)
        line2 = ogr.CreateGeometryFromWkt("LINESTRING (-80.5 40.5, -81.0 41.0)")
        feat2 = ogr.Feature(layer.GetLayerDefn())
        feat2.SetGeometry(line2)
        feat2.SetField("downstream_id", 0)
        layer.CreateFeature(feat2)
        ds = None

        # Mode 0: first feature last vertex
        lon0, lat0 = get_outlet_location(shp_path, iFlag_furthest_downstream=0)
        self.assertAlmostEqual(lon0, -80.5)
        self.assertAlmostEqual(lat0, 40.5)

        # Mode 1: terminal reach (downstream_id = 0) last vertex
        lon1, lat1 = get_outlet_location(shp_path, iFlag_furthest_downstream=1)
        self.assertAlmostEqual(lon1, -81.0)
        self.assertAlmostEqual(lat1, 41.0)


if __name__ == "__main__":
    unittest.main()
