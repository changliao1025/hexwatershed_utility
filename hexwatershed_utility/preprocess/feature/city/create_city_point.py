import os, sys
from osgeo import ogr, osr, gdal

from pyearth.system.define_global_variables import *
from pyearth.toolbox.reader.text_reader_string import text_reader_string

sWorkspace_data = '/compyfs/liao313/00raw/city/simplemaps_worldcities_basicv1.77'
sFilename_city = 'worldcities.csv'

sFilename = os.path.join(sWorkspace_data, sFilename_city)

aDatasets = text_reader_string(sFilename, iSkipline_in=1, cDelimiter_in=',')

aLongitude = aDatasets[:, 3]
aLatitude = aDatasets[:, 2]
aCity = aDatasets[:, 0]
aPopulation = aDatasets[:, 9]

nrow = len(aDatasets)

sFilename_out = os.path.join(sWorkspace_data, 'large_cities.geojson')
pDriver = ogr.GetDriverByName('GeoJSON')

#define spatial reference
srs = osr.SpatialReference()
srs.ImportFromEPSG(4326)

if os.path.exists(sFilename_out):
    pDriver.DeleteDataSource(sFilename_out)

pDataSource = pDriver.CreateDataSource(sFilename_out)
pLayer = pDataSource.CreateLayer('city', srs, ogr.wkbPoint)
#what attribute should be added?
pField = ogr.FieldDefn('city', ogr.OFTString)
pLayer.CreateField(pField)
pField = ogr.FieldDefn('latitude', ogr.OFTReal)
pLayer.CreateField(pField)
pField = ogr.FieldDefn('longitude', ogr.OFTReal)
pLayer.CreateField(pField)
pField = ogr.FieldDefn('population', ogr.OFTReal)
pLayer.CreateField(pField)
pLayerDefn = pLayer.GetLayerDefn()
pFeature = ogr.Feature(pLayerDefn)

for i in range(1, nrow):
    sName = aCity[i]
    dLatitude = float(aLatitude[i])
    dLongitude = float(aLongitude[i])

    sPopulation = aPopulation[i]
    if sPopulation == '':
        continue

    dPopulation = float(sPopulation)
    if dPopulation < 1.0E6:
        continue
    #create point geometry
    pPoint = ogr.Geometry(ogr.wkbPoint)
    pPoint.AddPoint(dLongitude, dLatitude)
    pFeature.SetGeometry(pPoint)
    pFeature.SetField('city', sName)
    pFeature.SetField('latitude', dLatitude)
    pFeature.SetField('longitude', dLongitude)
    pFeature.SetField('population', dPopulation)
    pLayer.CreateFeature(pFeature)
    pass

pDataSource.Destroy()

print('finished')

