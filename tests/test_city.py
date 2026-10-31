import os
import tempfile
import unittest

from osgeo import ogr, osr

from hexwatershed_utility.preprocess.feature.city import create_city_mask, create_city_point


class TestCity(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_file_not_found(self):
        with self.assertRaises(FileNotFoundError):
            create_city_point("non_existent.csv", "out.geojson")
        with self.assertRaises(FileNotFoundError):
            create_city_mask("non_existent.geojson", "out.tif")

    def test_create_city_point_from_csv(self):
        # Create sample CSV with header
        csv_path = os.path.join(self.temp_dir.name, "cities.csv")
        with open(csv_path, "w", encoding="utf-8") as f:
            f.write("city_name,country,latitude,longitude,population\n")
            f.write("MegaCity,CountryA,40.7128,-74.0060,8000000\n")
            f.write("SmallTown,CountryB,34.0522,-118.2437,50000\n")

        out_geojson = os.path.join(self.temp_dir.name, "cities.geojson")
        res = create_city_point(
            csv_path,
            out_geojson,
            dPopulation_min=1.0e6,
            sCol_city="city_name",
            sCol_lat="latitude",
            sCol_lon="longitude",
            sCol_pop="population",
        )
        self.assertTrue(os.path.exists(res))

        ds = ogr.Open(out_geojson)
        layer = ds.GetLayer(0)
        self.assertEqual(layer.GetFeatureCount(), 1)
        feat = layer.GetFeature(0)
        self.assertEqual(feat.GetField("city"), "MegaCity")
        self.assertAlmostEqual(feat.GetField("population"), 8000000.0)
        ds = None

    def test_create_city_point_from_vector(self):
        # Create sample vector dataset
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)

        driver = ogr.GetDriverByName("GeoJSON")
        vector_path = os.path.join(self.temp_dir.name, "raw_cities.geojson")
        ds = driver.CreateDataSource(vector_path)
        layer = ds.CreateLayer("cities", srs, ogr.wkbPoint)
        layer.CreateField(ogr.FieldDefn("NAME", ogr.OFTString))
        layer.CreateField(ogr.FieldDefn("POP_MAX", ogr.OFTReal))

        # Big city
        pt1 = ogr.Geometry(ogr.wkbPoint)
        pt1.AddPoint(116.4, 39.9)
        f1 = ogr.Feature(layer.GetLayerDefn())
        f1.SetGeometry(pt1)
        f1.SetField("NAME", "CapitalCity")
        f1.SetField("POP_MAX", 20000000.0)
        layer.CreateFeature(f1)

        # Small town
        pt2 = ogr.Geometry(ogr.wkbPoint)
        pt2.AddPoint(116.5, 40.0)
        f2 = ogr.Feature(layer.GetLayerDefn())
        f2.SetGeometry(pt2)
        f2.SetField("NAME", "Village")
        f2.SetField("POP_MAX", 5000.0)
        layer.CreateFeature(f2)
        ds = None

        out_geojson = os.path.join(self.temp_dir.name, "filtered_cities.geojson")
        res = create_city_point(vector_path, out_geojson, dPopulation_min=1.0e6)
        self.assertTrue(os.path.exists(res))

        ds_out = ogr.Open(out_geojson)
        layer_out = ds_out.GetLayer(0)
        self.assertEqual(layer_out.GetFeatureCount(), 1)
        feat_out = layer_out.GetFeature(0)
        self.assertEqual(feat_out.GetField("city"), "CapitalCity")
        ds_out = None


if __name__ == "__main__":
    unittest.main()
