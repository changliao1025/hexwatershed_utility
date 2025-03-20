import os
from osgeo import ogr, gdal, osr
sFolder  = '/qfs/people/liao313/data/hexwatershed/greatlakes/vector/hydrology'

aFilename = list()
aFilename.append( 'lake_erie.geojson' )
aFilename.append( 'lake_huron.geojson' )
aFilename.append( 'lake_michigan.geojson' )
aFilename.append( 'lake_ontario.geojson' )
aFilename.append( 'lake_superior.geojson' )

pDriver_geojson = ogr.GetDriverByName('GeoJSON')
for sFilename in aFilename:
    sFilename = os.path.join( sFolder, sFilename )

    if os.path.exists( sFilename ) == False:
        print( sFilename )

    pDatasetIn = pDriver_geojson.Open(sFilename, 0)
    pLayerIn = pDatasetIn.GetLayer()
    pSpatialRef = pLayerIn.GetSpatialRef()

    sFilename_new = sFilename.replace('.geojson', '_new.geojson')
    if os.path.exists( sFilename_new ):
        os.remove( sFilename_new )
    else:
        print(sFilename_new)

    pDatasetOut = pDriver_geojson.CreateDataSource(sFilename_new)

    # Create new pLayerIn with the same spatial reference as the original
    pLayerOut = pDatasetOut.CreateLayer('lake', pSpatialRef, ogr.wkbPolygon)
    pLayerOut.CreateField(ogr.FieldDefn('lakeid', ogr.OFTInteger)) #long type for high resolution

    pFeatureOut = ogr.Feature(pLayerOut.GetLayerDefn())
    #read all the features in the pLayerIn
    nFeature = pLayerIn.GetFeatureCount()

    pGeometry_union  = ogr.Geometry(ogr.wkbPolygon)
    for i in range(nFeature):
        pFeature = pLayerIn.GetFeature(i)
        #get pGeometry
        pGeometry = pFeature.GetGeometryRef()
        #get pGeometry type
        sGeometryName = pGeometry.GetGeometryName()
        if sGeometryName == 'POLYGON':
            ring = pGeometry.GetGeometryRef(0)
            sGeometryName1 = ring.GetGeometryName()
            if sGeometryName1 == 'POLYGON':
                ring1 = ring.GetGeometryRef(0)
                exterior  = ogr.Geometry(ogr.wkbPolygon)
                exterior.AddGeometry(ring1)
                pGeometry_union = pGeometry_union.Union(exterior)

        else:
            if sGeometryName == 'MULTIPOLYGON':
                ring = pGeometry.GetGeometryRef(0)
                sGeometryName1 = ring.GetGeometryName()
                if sGeometryName1 == 'POLYGON':
                    ring1 = ring.GetGeometryRef(0)
                    exterior  = ogr.Geometry(ogr.wkbPolygon)
                    exterior.AddGeometry(ring1)
                    pGeometry_union = pGeometry_union.Union(exterior)



    pFeatureOut.SetGeometry(pGeometry_union)
    pFeatureOut.SetField("lakeid", 1)
    pLayerOut.CreateFeature(pFeatureOut)
    pFeatureOut = None
    pDatasetOut = None #create a new pLayerIn






