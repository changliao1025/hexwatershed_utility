import os, sys
import numpy as np
from osgeo import ogr, osr, gdal
from tinyr import RTree

from pyearth.system.define_global_variables import *
from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates
from pyearth.gis.geometry.calculate_polygon_area import calculate_polygon_area
from pyflowline.classes.vertex import pyvertex
from pyflowline.formats.read_flowline import read_flowline_geojson
sWorkspace_data = '/compyfs/liao313/00raw/dam/GOODD_data/Data'

sFilename_in = os.path.join(sWorkspace_data, 'GOOD2_dams.geojson')
sFilename_out = os.path.join(sWorkspace_data, 'GOOD2_dams_snap.geojson')

sFilename_river= '/qfs/people/liao313/data/hexwatershed/global/vector/river_networks_wo_lakes.geojson'
if not os.path.exists(sFilename_river):
    print('river network does not exist')
    sys.exit()

pDriver = ogr.GetDriverByName('GeoJSON')

pDataSource_dam = pDriver.Open(sFilename_in, 0)
pLayer_dam = pDataSource_dam.GetLayer()
nDam = pLayer_dam.GetFeatureCount()

pSpatialReference = osr.SpatialReference()
pSpatialReference.ImportFromEPSG(4326)

if os.path.exists(sFilename_out):
    pDriver.DeleteDataSource(sFilename_out)

pDataSource_out = pDriver.CreateDataSource(sFilename_out)
pLayer_out = pDataSource_out.CreateLayer('dam', pSpatialReference, ogr.wkbPoint)
pLayer_out.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
pLayer_out.CreateField(ogr.FieldDefn('dLongitude_degree', ogr.OFTReal))
pLayer_out.CreateField(ogr.FieldDefn('dLatitude_degree', ogr.OFTReal))
pLayer_out.CreateField(ogr.FieldDefn('drainage', ogr.OFTReal))
pLayerDefn = pLayer_out.GetLayerDefn()
pFeature_out = ogr.Feature(pLayerDefn)

pDataSource_river = pDriver.Open(sFilename_river, 0)
pLayer_river = pDataSource_river.GetLayer()
nFlowline = pLayer_river.GetFeatureCount()

aFlowline, pProjection_geojson = read_flowline_geojson( sFilename_river )

index_flowline = RTree(max_cap=5, min_cap=2)
for i in range(nFlowline):
    pBound= aFlowline[i].pBound
    index_flowline.insert(i, pBound)  #
    pass

#read dam location
iDam_id = 1
dBuffer = 0.01 #degree, approximately 1 km
for i in range(nDam):
    pFeature_dam = pLayer_dam.GetFeature(i)
    point = pFeature_dam.GetGeometryRef()
    dDrainage = pFeature_dam.GetField('drainage')
    #get point coordinates using gdal api
    x = point.GetX()
    y = point.GetY()
    point0= dict()
    point0['dLongitude_degree'] = x
    point0['dLatitude_degree'] = y
    pVertex=pyvertex(point0)
    #aIntersect = list(index_flowline.search_surrounding([x,y]))
    left =   x - dBuffer
    right =  x + dBuffer
    bottom = y - dBuffer
    top =    y + dBuffer
    pBound= (left, bottom, right, top)
    aIntersect = list(index_flowline.search( pBound )  )
    if len(aIntersect) == 0:
        #print('dam is not close to any river')
        pPoint = ogr.Geometry(ogr.wkbPoint)
        pPoint.AddPoint(x, y)
        #create a new feature
        pFeature_out.SetGeometry(pPoint)
        pFeature_out.SetField('id', iDam_id)
        pFeature_out.SetField('dLongitude_degree', x)
        pFeature_out.SetField('dLatitude_degree', y)
        pFeature_out.SetField('drainage', dDrainage)
        pLayer_out.CreateFeature(pFeature_out)
        iDam_id = iDam_id + 1
        pass
    else:
        dDistance_min = 1.0E8
        lIndex_closest = -1
        if i == 48:
            pass
        for k in aIntersect:
            #calculate the distance from point to the flowline
            pFlowline = aFlowline[k]
            distance, pVertex_out = pFlowline.calculate_distance_to_vertex(pVertex)
            if distance < dDistance_min:
                dDistance_min = distance
                lIndex_closest = k
                pVertex_closest = pVertex_out
                pass
            pass

        if dDistance_min < 3.0E3:
            #save the dam location
            #create point geometry
            pPoint = ogr.Geometry(ogr.wkbPoint)
            pPoint.AddPoint(pVertex_closest.dLongitude_degree, pVertex_closest.dLatitude_degree)
            #create a new feature
            pFeature_out.SetGeometry(pPoint)
            pFeature_out.SetField('id', iDam_id)
            pFeature_out.SetField('dLongitude_degree', pVertex_closest.dLongitude_degree)
            pFeature_out.SetField('dLatitude_degree', pVertex_closest.dLatitude_degree)
            pFeature_out.SetField('drainage', dDrainage)
            pLayer_out.CreateFeature(pFeature_out)
            iDam_id = iDam_id + 1
            pass
        else:
            print('dam is too far from any river')
            print(i, x, y, 'distance is ', dDistance_min)

            pPoint = ogr.Geometry(ogr.wkbPoint)
            pPoint.AddPoint(x, y)
            #create a new feature
            pFeature_out.SetGeometry(pPoint)
            pFeature_out.SetField('id', iDam_id)
            pFeature_out.SetField('dLongitude_degree', x)
            pFeature_out.SetField('dLatitude_degree', y)
            pFeature_out.SetField('drainage', dDrainage)
            pLayer_out.CreateFeature(pFeature_out)
            iDam_id = iDam_id + 1
            pass



    pass

pDataSource_dam.Destroy()
pDataSource_river.Destroy()
pDataSource_out.Destroy()
print('finished!')
