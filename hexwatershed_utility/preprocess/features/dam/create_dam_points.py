import os, sys
import numpy as np
from osgeo import ogr, osr, gdal

from pyearth.system.define_global_variables import *
from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates
from pyearth.gis.geometry.calculate_polygon_area import calculate_polygon_area

sWorkspace_data = '/compyfs/liao313/00raw/dam/GOODD_data/Data'

sFilename_dam = os.path.join(sWorkspace_data, 'GOOD2_dams.shp')

sFilename_catchment = os.path.join(sWorkspace_data, 'GOOD2_catchments.shp')

sFilename_out = os.path.join(sWorkspace_data, 'GOOD2_dams.geojson')

pDriver_shp = ogr.GetDriverByName('ESRI Shapefile')
pDriver_json = ogr.GetDriverByName('GeoJSON')


#define spatial reference
srs = osr.SpatialReference()
srs.ImportFromEPSG(4326)
pProjection = srs.ExportToWkt()

if os.path.exists(sFilename_out):
    pDriver_json.DeleteDataSource(sFilename_out)

pDataSource = pDriver_json.CreateDataSource(sFilename_out)
pLayer = pDataSource.CreateLayer('dam', srs, ogr.wkbPoint)
#what attribute should be added?
pField = ogr.FieldDefn('dam', ogr.OFTString)
pLayer.CreateField(pField)
pField = ogr.FieldDefn('latitude', ogr.OFTReal)
pLayer.CreateField(pField)
pField = ogr.FieldDefn('longitude', ogr.OFTReal)
pLayer.CreateField(pField)
pField = ogr.FieldDefn('drainage', ogr.OFTReal)
pLayer.CreateField(pField)
pLayerDefn = pLayer.GetLayerDefn()
pFeature = ogr.Feature(pLayerDefn)

#open the shapefile
pDataSource_dam = pDriver_shp.Open(sFilename_dam, 0)
pLayer_dam = pDataSource_dam.GetLayer()
#open the shapefile
pDataSource_catchment = pDriver_shp.Open(sFilename_catchment, 0)
pLayer_catchment = pDataSource_catchment.GetLayer()

#get the spatial refernce of the shapefile
pSpatialRef_dam = pLayer_dam.GetSpatialRef()
pSpatialRef_catchment = pLayer_catchment.GetSpatialRef()

#convert to wkt projection
srs_wkt = pSpatialRef_dam.ExportToWkt()
srs_wkt_catchment = pSpatialRef_catchment.ExportToWkt()
if srs_wkt != pProjection:
    iProjection_dam = 1
else:
    iProjection_dam = 0

if srs_wkt_catchment != pProjection:
    iProjection_catchment = 1
else:
    iProjection_catchment = 0


#read all the features in the dam shapefile
nFeature_dam = pLayer_dam.GetFeatureCount()
nFeature_catchment = pLayer_catchment.GetFeatureCount()

aLongitude = []
aLatitude = []
aDam_id = []
for i in range(nFeature_dam):
    pFeature_dam = pLayer_dam.GetFeature(i)
    pGeometry_dam = pFeature_dam.GetGeometryRef()
    if iProjection_dam == 1:
        pGeometry_dam.TransformTo(srs)

    # List all attributes of the dam feature
    #pDefn = pFeature_dam.GetDefnRef()
    #nFields = pDefn.GetFieldCount()
    #print(f"Feature {i} attributes:")
    #for j in range(nFields):
    #    fieldDefn = pDefn.GetFieldDefn(j)
    #    fieldName = fieldDefn.GetNameRef()
    #    fieldValue = pFeature_dam.GetField(j)
    #    print(f"  {fieldName}: {fieldValue}")

    pPoint_dam = pGeometry_dam.GetPoint()
    dLongitude = pPoint_dam[0]
    dLatitude = pPoint_dam[1]
    #create point geometry
    #get the dam id
    iDam_id = pFeature_dam.GetField('DAM_ID')
    aLongitude.append(dLongitude)
    aLatitude.append(dLatitude)
    aDam_id.append(iDam_id)

aDam_id = np.array(aDam_id)
#read all the features in the catchment shapefile
for i in range(nFeature_catchment):
    pFeature_catchment = pLayer_catchment.GetFeature(i)
    pGeometry_catchment = pFeature_catchment.GetGeometryRef()
    if iProjection_catchment == 1:
        pGeometry_catchment.TransformTo(srs)
    #create point geometry
    #get the dam id
    iDam_id = pFeature_catchment.GetField('DAM_ID')
    #use this id to find the corresponding dam
    dummy= np.where(aDam_id == iDam_id)
    if len(dummy[0]) == 0:
        continue
    else:
        iIndex = dummy[0][0]

    dLongitude = aLongitude[iIndex]
    dLatitude = aLatitude[iIndex]

    aCoords_gcs = get_geometry_coordinates(pGeometry_catchment)
    sGeometry_type = pGeometry_catchment.GetGeometryName()
    if sGeometry_type =='MULTIPOLYGON':
        #pick the one with larger area in the list of coords
        aArea = []
        for aCoords in aCoords_gcs:
            dArea = calculate_polygon_area(aCoords[:,0], aCoords[:,1])
            aArea.append(dArea)
        dArea = max(aArea)
    else:
        dArea = calculate_polygon_area(aCoords_gcs[:,0], aCoords_gcs[:,1])

    dArea = dArea/1.0E6 #km2
    if dArea < 1.0E4:
        continue

    #create point geometry
    pPoint = ogr.Geometry(ogr.wkbPoint)
    pPoint.AddPoint(dLongitude, dLatitude)
    #create a new feature
    pFeature.SetGeometry(pPoint)
    pFeature.SetField('dam', iDam_id)
    pFeature.SetField('latitude', dLatitude)
    pFeature.SetField('longitude', dLongitude)
    pFeature.SetField('drainage', dArea)
    pLayer.CreateFeature(pFeature)
    pass

#close the files
pDataSource.Destroy()
pDataSource_dam.Destroy()
pDataSource_catchment.Destroy()

print('finished')
