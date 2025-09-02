import os
from osgeo import ogr, gdal, osr
from pyearth.toolbox.data.geoparquet.convert_geojson_to_geoparquet import convert_geojson_to_geoparquet

sFolder  = '/qfs/people/liao313/data/hexwatershed/greatlakes/vector/hydrology'
pDriver_geojson = ogr.GetDriverByName('GeoJSON')
#read the whole domain
sFilename = os.path.join( sFolder, 'basins_boundary.geojson' )
pDataset_whole = pDriver_geojson.Open(sFilename, 0)
pLayer_whole = pDataset_whole.GetLayer()
pSpatialRef = pLayer_whole.GetSpatialRef()
nFeature_whole = pLayer_whole.GetFeatureCount()

sFilename_out = os.path.join( sFolder, 'lake_difference.geojson' )
if os.path.exists( sFilename_out ):
    os.remove( sFilename_out )

pDatesetOut = pDriver_geojson.CreateDataSource(sFilename_out)
pLayerOut = pDatesetOut.CreateLayer('lake', pSpatialRef, ogr.wkbPolygon)
pLayerOut.CreateField(ogr.FieldDefn('polygonid', ogr.OFTInteger)) #l
pFeatureOut = ogr.Feature(pLayerOut.GetLayerDefn())

pFeature = pLayer_whole.GetFeature(0)
pGeometry = pFeature.GetGeometryRef()
sGeometryName = pGeometry.GetGeometryName()
pPolygonOut = ogr.Geometry(ogr.wkbPolygon)

if sGeometryName == 'MULTIPOLYGON':
    ring = pGeometry.GetGeometryRef(0)
    sGeometryName1 = ring.GetGeometryName()
    if sGeometryName1 == 'POLYGON':
        ring1 = ring.GetGeometryRef(0)
        pPolygonOut.AddGeometry(ring1)

dBuffer_threshold = 0.001

aFilename = list()
aFilename.append( 'lake_erie_new.geojson' )
aFilename.append( 'lake_huron_new.geojson' )
aFilename.append( 'lake_michigan_new.geojson' )
aFilename.append( 'lake_ontario_new.geojson' )
aFilename.append( 'lake_superior_new.geojson' )

for sFilename in aFilename:
    sFilename = os.path.join( sFolder, sFilename )

    if os.path.exists( sFilename ) == False:
        print( sFilename )

    pDataset_base = pDriver_geojson.Open(sFilename, 0)
    pLayerIn = pDataset_base.GetLayer()
    pSpatialRef = pLayerIn.GetSpatialRef()
    nFeature = pLayerIn.GetFeatureCount()
    for i in range(nFeature):
        pFeature = pLayerIn.GetFeature(i)
        #get pGeometry
        pGeometry = pFeature.GetGeometryRef()
        #get pGeometry type
        sGeometryName = pGeometry.GetGeometryName()
        if sGeometryName == 'POLYGON':
            pPolygon = pGeometry.GetGeometryRef(0)
            sGeometryName = pPolygon.GetGeometryName()
            if sGeometryName == 'LINEARRING':
                pPolygonOut.AddGeometry(pPolygon)

            #create a hole
        else:
            if sGeometryName == 'MULTIPOLYGON':
                iCount = pGeometry.GetGeometryCount()
                for j in range(iCount):
                    pPolygon = pGeometry.GetGeometryRef(j)
                    sGeometryName1 = pPolygon.GetGeometryName()
                    if sGeometryName1 == 'POLYGON':
                        ring_inner = pPolygon.GetGeometryRef(0)
                        sGeometryName2 = ring_inner.GetGeometryName()
                        if sGeometryName2 == 'LINEARRING':
                            pPolygonOut.AddGeometry(ring_inner)

            print('sGeometryName:', sGeometryName)

pFeatureOut.SetGeometry(pPolygonOut)
pFeatureOut.SetField('polygonid', 1)
pLayerOut.CreateFeature(pFeatureOut)

#flush the output

pLayerOut = None
pDatesetOut = None
pDataset_new = None

#convert to partquet
sFilename_out2 = sFilename_out.replace('.geojson', '.parquet')
convert_geojson_to_geoparquet(sFilename_out, sFilename_out2)



