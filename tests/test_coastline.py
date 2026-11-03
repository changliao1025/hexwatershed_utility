import os
import tempfile
import unittest

from osgeo import ogr, osr
from pyearth.toolbox.management.vector.remove_internal_polygon import remove_internal_polygon

from hexwatershed_utility.preprocess.feature.coastline import (
    create_land_ocean_mask_from_hydrobasins,
    create_land_ocean_mask_from_naturalearth,
    extend_land_mask_to_rivers,
    geometries_bbox_overlap,
)


class TestCoastline(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_callables_exist(self):
        self.assertTrue(callable(create_land_ocean_mask_from_hydrobasins))
        self.assertTrue(callable(create_land_ocean_mask_from_naturalearth))
        self.assertTrue(callable(extend_land_mask_to_rivers))
        self.assertTrue(callable(geometries_bbox_overlap))

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
            create_land_ocean_mask_from_hydrobasins("/tmp/out")

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

    def test_simplify_vector_geometry(self):
        from hexwatershed_utility.preprocess.feature.coastline import simplify_vector_geometry

        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)

        driver = ogr.GetDriverByName("GeoJSON")
        in_geojson = os.path.join(self.temp_dir.name, "complex_poly.geojson")
        ds = driver.CreateDataSource(in_geojson)
        layer = ds.CreateLayer("land", srs, ogr.wkbPolygon)

        # Polygon with many intermediate redundant collinear points along an edge
        wkt_complex = "POLYGON ((0 0, 1 0, 2 0, 3 0, 4 0, 5 0, 5 5, 0 5, 0 0))"
        geom = ogr.CreateGeometryFromWkt(wkt_complex)
        feat = ogr.Feature(layer.GetLayerDefn())
        feat.SetGeometry(geom)
        layer.CreateFeature(feat)
        ds = None

        out_geojson = os.path.join(self.temp_dir.name, "simplified_poly.geojson")
        simplify_vector_geometry(in_geojson, out_geojson, dTolerance=0.01)
        self.assertTrue(os.path.exists(out_geojson))

        ds_out = ogr.Open(out_geojson)
        layer_out = ds_out.GetLayer(0)
        feat_out = layer_out.GetFeature(0)
        ring_out = feat_out.GetGeometryRef().GetGeometryRef(0)
        # Originally 9 points, should simplify down to 5 points (a rectangle)
        self.assertLess(ring_out.GetPointCount(), 9)
        self.assertEqual(ring_out.GetPointCount(), 5)
        ds_out = None

    def test_create_land_ocean_mask_from_naturalearth_vector_only(self):
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)

        driver = ogr.GetDriverByName("GeoJSON")
        in_geojson = os.path.join(self.temp_dir.name, "test_input_land.geojson")
        ds = driver.CreateDataSource(in_geojson)
        layer = ds.CreateLayer("land", srs, ogr.wkbPolygon)

        # Polygon of reasonable size
        wkt = "POLYGON ((10 10, 12 10, 12 12, 10 12, 10 10))"
        geom = ogr.CreateGeometryFromWkt(wkt)
        feat = ogr.Feature(layer.GetLayerDefn())
        feat.SetGeometry(geom)
        layer.CreateFeature(feat)
        ds = None

        sFilename_buf_out = os.path.join(self.temp_dir.name, "test_buffered.geojson")
        tif_out, vec_out = create_land_ocean_mask_from_naturalearth(
            sWorkspace_output=self.temp_dir.name,
            sFilename_land_vector_in=in_geojson,
            dThreshold_area_island=0.0,
            iFlag_simplify_coastline=True,
            dDistance_buffer_meter=5000.0,
            sFilename_vector_out=sFilename_buf_out,
            iFlag_convert_raster=False,
        )

        self.assertIsNone(tif_out)
        self.assertEqual(vec_out, sFilename_buf_out)
        self.assertTrue(os.path.exists(sFilename_buf_out))
        ds_buf = ogr.Open(sFilename_buf_out)
        self.assertGreater(ds_buf.GetLayer(0).GetFeatureCount(), 0)
        ds_buf = None


if __name__ == "__main__":
    unittest.main()
