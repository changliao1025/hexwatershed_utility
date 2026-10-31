import os
import tempfile
import unittest

from osgeo import ogr, osr

from hexwatershed_utility.preprocess.feature.lake import create_lake_mask, extract_lake_polygon


class TestLake(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_extract_lake_polygon_file_not_found(self):
        with self.assertRaises(FileNotFoundError):
            extract_lake_polygon("non_existent.shp", "out.geojson")

    def test_create_lake_mask_file_not_found(self):
        with self.assertRaises(FileNotFoundError):
            create_lake_mask("non_existent.geojson", "out.tif")

    def test_extract_lake_polygon_basic(self):
        # Create a sample GeoJSON with two polygons: one large and one small
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)

        driver = ogr.GetDriverByName("GeoJSON")
        shp_path = os.path.join(self.temp_dir.name, "sample_lakes.geojson")
        ds = driver.CreateDataSource(shp_path)
        layer = ds.CreateLayer("lakes", srs, ogr.wkbPolygon)

        layer.CreateField(ogr.FieldDefn("name", ogr.OFTString))
        layer.CreateField(ogr.FieldDefn("lake_type", ogr.OFTInteger))

        # Big lake (~10 degrees x 10 degrees, > 5000 km2)
        poly_big = ogr.CreateGeometryFromWkt("POLYGON ((-85 45, -75 45, -75 50, -85 50, -85 45))")
        feat = ogr.Feature(layer.GetLayerDefn())
        feat.SetGeometry(poly_big)
        feat.SetField("name", "Big Lake")
        feat.SetField("lake_type", 1)
        layer.CreateFeature(feat)

        # Tiny lake (< 5000 km2)
        poly_small = ogr.CreateGeometryFromWkt(
            "POLYGON ((-80.1 45.1, -80.0 45.1, -80.0 45.2, -80.1 45.2, -80.1 45.1))"
        )
        feat = ogr.Feature(layer.GetLayerDefn())
        feat.SetGeometry(poly_small)
        feat.SetField("name", "Tiny Lake")
        feat.SetField("lake_type", 2)
        layer.CreateFeature(feat)

        ds = None

        out_geojson = os.path.join(self.temp_dir.name, "extracted_lakes.geojson")
        res = extract_lake_polygon(shp_path, out_geojson, dThreshold_area=5000.0)
        self.assertTrue(os.path.exists(res))

        # Check that only the big lake was extracted
        ds_out = ogr.Open(out_geojson)
        layer_out = ds_out.GetLayer(0)
        self.assertEqual(layer_out.GetFeatureCount(), 1)
        feat_out = layer_out.GetFeature(0)
        self.assertGreater(feat_out.GetField("area"), 5000.0)
        ds_out = None


if __name__ == "__main__":
    unittest.main()
