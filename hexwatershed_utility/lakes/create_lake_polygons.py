import os, sys
import numpy as np
from osgeo import ogr, osr, gdal

from pyearth.system.define_global_variables import *
from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates
from pyearth.gis.geometry.calculate_polygon_area import calculate_polygon_area

sWorkspace_data = '/compyfs/liao313/00raw/hydrology/hydroshed/hydrolake/HydroLAKES_polys_v10_shp'
sWorkspace_data_out  = '/compyfs/liao313/00raw/hydrology/hydroshed/hydrolake/'



sFilename_out = os.path.join(sWorkspace_data_out, 'global_lakes.geojson')

pDriver_shp = ogr.GetDriverByName('ESRI Shapefile')
pDriver_json = ogr.GetDriverByName('GeoJSON')

sFilename_in = os.path.join(sWorkspace_data, 'HydroLAKES_polys_v10.shp')

#open the shapefile
pDataSource_in = pDriver_shp.Open(sFilename_in, 0)
pLayer_in = pDataSource_in.GetLayer()
#get spatial reference
pSpatialRef_in = pLayer_in.GetSpatialRef()
srs_wkt = pSpatialRef_in.ExportToWkt()

if os.path.exists(sFilename_out):
    pDriver_json.DeleteDataSource(sFilename_out)

pDataSource_out = pDriver_json.CreateDataSource(sFilename_out)
srs = osr.SpatialReference()
srs.ImportFromEPSG(4326)
pLayer_out = pDataSource_out.CreateLayer('lake', srs, ogr.wkbPolygon)
pLayerDefn = pLayer_out.GetLayerDefn()
#add area field
pField = ogr.FieldDefn('area', ogr.OFTReal)
pLayer_out.CreateField(pField)
pFeature_out = ogr.Feature(pLayerDefn)

sProjection = srs.ExportToWkt()

if srs_wkt != sProjection:
    iProjection = 1
else:
    iProjection = 0


# Set attribute filter to only process large lakes

dThreshold_area = 5000.0 #km2

#convert to string using scientific notation
sThreshold_area = "{:.1E}".format(dThreshold_area)

sFilter = "Lake_area > " + sThreshold_area

pLayer_in.SetAttributeFilter(sFilter) #default unit is km2

# Get the number of filtered features
filtered_feature_count = pLayer_in.GetFeatureCount()
print(f"Number of filtered features: {filtered_feature_count}")

for pFeature_in in pLayer_in:
    pGeometry_in = pFeature_in.GetGeometryRef()

    if iProjection == 1:
        pGeometry_in.TransformTo(srs)

    #pGeometry_out = pGeometry_in.Clone()
    exterior_ring = pGeometry_in.GetGeometryRef(0)
    pGeometry_out = ogr.Geometry(ogr.wkbPolygon)
    pGeometry_out.AddGeometry(exterior_ring)

    aCoords_gcs = get_geometry_coordinates(pGeometry_in)
    sGeometry_type = pGeometry_in.GetGeometryName()

    dArea = calculate_polygon_area(aCoords_gcs[:,0], aCoords_gcs[:,1])

    dArea = dArea/1.0e6 #convert to km2
    if dArea < dThreshold_area: #we only keep large lake that has area larger than 1 km2
        continue

    pFeature_out.SetGeometry(pGeometry_out)
    pFeature_out.SetField('area', dArea)
    pLayer_out.CreateFeature(pFeature_out)

pDataSource_out.Destroy()
pDataSource_in.Destroy()
print('finished')
