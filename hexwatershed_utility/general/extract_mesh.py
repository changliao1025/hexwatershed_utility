import os
from osgeo import ogr, gdal, osr
from pyearth.toolbox.data.geoparquet.convert_geojson_to_geoparquet import convert_geojson_to_geoparquet

sFilename_mesh  = '/compyfs/liao313/04model/pyflowline/greatlakes/pyflowline20230701005/mpas.geojson'

pDriver_geojson = ogr.GetDriverByName('GeoJSON')
#read the whole domain
pDataset_mesh = pDriver_geojson.Open(sFilename_mesh, 0)
pLayer_mesh = pDataset_mesh.GetLayer()
pSpatialRef = pLayer_mesh.GetSpatialRef()
nFeature_mesh = pLayer_mesh.GetFeatureCount()


sFilename_out = '/compyfs/liao313/04model/pyflowline/greatlakes/pyflowline20230701005/mpas_mesh_extract.geojson'


if os.path.exists( sFilename_out ):
    os.remove( sFilename_out )

pDatesetOut = pDriver_geojson.CreateDataSource(sFilename_out)
pLayerOut = pDatesetOut.CreateLayer('cell', pSpatialRef, ogr.wkbPolygon)
pLayerOut.CreateField(ogr.FieldDefn('cellid', ogr.OFTInteger)) #l
pFeatureOut = ogr.Feature(pLayerOut.GetLayerDefn())

#read the lake boundary



sFolder  = '/qfs/people/liao313/data/hexwatershed/greatlakes/vector/hydrology'
aFilename = list()
aFilename.append( 'lake_erie_new.geojson' )
aFilename.append( 'lake_huron_new.geojson' )
aFilename.append( 'lake_michigan_new.geojson' )
aFilename.append( 'lake_ontario_new.geojson' )
aFilename.append( 'lake_superior_new.geojson' )
nFile = len(aFilename)
aGeometry_lake = list()
for sFilename in aFilename:
    sFilename = os.path.join( sFolder, sFilename )
    if os.path.exists( sFilename ) == False:
        print( sFilename )

    pGeometry_lake = ogr.Geometry(ogr.wkbPolygon)

    pDataset_base = pDriver_geojson.Open(sFilename, 0)
    pLayerIn = pDataset_base.GetLayer()
    pSpatialRef = pLayerIn.GetSpatialRef()
    nFeature = pLayerIn.GetFeatureCount()
    for i in range(nFeature):
        pFeature = pLayerIn.GetFeature(i)
        pGeometry = pFeature.GetGeometryRef()
        sGeometryName = pGeometry.GetGeometryName()
        if sGeometryName == 'POLYGON':
            pPolygon = pGeometry.GetGeometryRef(0)
            ring = ogr.Geometry(ogr.wkbLinearRing)
            nVertex = pPolygon.GetPointCount()
            for k in range(nVertex):
                x1 = pPolygon.GetX(k)
                y1 = pPolygon.GetY(k)
                ring.AddPoint(x1, y1)
            pGeometry_lake.AddGeometry(ring)
            aGeometry_lake.append(pGeometry_lake)
        else:
            if sGeometryName == 'MULTIPOLYGON':
                iCount = pGeometry.GetGeometryCount()
                for j in range(iCount):
                    pPolygon = pGeometry.GetGeometryRef(j)
                    sGeometryName1 = pPolygon.GetGeometryName()
                    if sGeometryName1 == 'POLYGON':
                        pPolygon1 = pPolygon.GetGeometryRef(0)
                        ring = ogr.Geometry(ogr.wkbLinearRing)
                        nVertex = pPolygon1.GetPointCount()
                        for k in range(nVertex):
                            x1 = pPolygon1.GetX(k)
                            y1 = pPolygon1.GetY(k)
                            ring.AddPoint(x1, y1)
                        pGeometry_lake.AddGeometry(ring)
                        aGeometry_lake.append(pGeometry_lake)

cellid = 1
for i in range(nFeature_mesh):
    pFeature_mesh = pLayer_mesh.GetFeature(i)
    pGeometry_mesh = pFeature_mesh.GetGeometryRef()
    sGeometryName_mesh = pGeometry_mesh.GetGeometryName()
    if sGeometryName_mesh == 'POLYGON':
        iFlag_within = False
        for j in range(nFile):
            pGeometry_lake = aGeometry_lake[j]
            iFlag_within = pGeometry_mesh.Within(pGeometry_lake)
            #iFlag_within = pGeometry_mesh.Intersects( pGeometry_lake )
            if iFlag_within == True:
                break

        if iFlag_within == False:
            pFeatureOut.SetGeometry(pGeometry_mesh)
            pFeatureOut.SetField('cellid', cellid)
            pLayerOut.CreateFeature(pFeatureOut)
            cellid = cellid + 1

pFeatureOut = None
pDatesetOut = None
print('done')



