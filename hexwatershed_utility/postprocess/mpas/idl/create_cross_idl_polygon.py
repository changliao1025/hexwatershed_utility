import os
import numpy as np
from osgeo import ogr, gdal, osr
from shapely.geometry import Polygon
from shapely.ops import transform
from shapely.wkt import loads
import pyproj
from pyearth.system.define_global_variables import *
from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates

from pyearth.gis.geometry.calculate_intersect_on_great_circle import find_great_circle_intersection

def split_polygon_at_idl_manual(aCoord_gcs):

    #find out the two index where the edge crosses the IDL
    nPoint = len(aCoord_gcs)

    aIndex = []

    #this algo require the coordinates to be in CCW order

    for i in range(nPoint - 1):
        dLongitude = aCoord_gcs[i,0]
        dLongitude_next = aCoord_gcs[i + 1,0]
        if dLongitude > 0 and  dLongitude < 180.0 and dLongitude_next < 0:
            aIndex.append(i)
            continue
        if dLongitude < 0 and dLongitude_next > 0:
            aIndex.append(i)
            continue

    if len(aIndex) != 2:
        print('Warning: no intersection found')
        return

    #get the two intersection points
    lon1= aCoord_gcs[aIndex[0],0]
    lat1= aCoord_gcs[aIndex[0],1]
    lon2= aCoord_gcs[aIndex[0] + 1,0]
    lat2= aCoord_gcs[aIndex[0] + 1,1]
    target_lon = 180.0
    d, dLat0 = find_great_circle_intersection(lon1, lat1, lon2, lat2, target_lon)

    lon1= aCoord_gcs[aIndex[1],0]
    lat1= aCoord_gcs[aIndex[1],1]
    if aIndex[1] == nPoint - 1:
        lon2= aCoord_gcs[0,0]
        lat2= aCoord_gcs[0,1]
    else:
        lon2= aCoord_gcs[aIndex[1] + 1,0]
        lat2= aCoord_gcs[aIndex[1] + 1,1]
    target_lon = 180.0
    d, dLat1 = find_great_circle_intersection(lon1, lat1, lon2, lat2, target_lon)

    #compare the two intersection points which is top and bottom
    if dLat0 > dLat1:
        dLat_top = dLat0
        dLat_bottom = dLat1
    else:
        dLat_top = dLat1
        dLat_bottom = dLat0

    aCoord_gcs_left = list()
    aCoord_gcs_right = list()

    #left part dLongitude > 0
    aIndex_dummy = np.array(aIndex)
    iFlag_added_right = 0
    iFlag_added_left = 0
    for i in range(nPoint):
        dLongitude = aCoord_gcs[i,0]
        if dLongitude > 0:
            if i <= np.min(aIndex_dummy):
                aCoord_gcs_left.append(aCoord_gcs[i])
                if i in aIndex:
                    if iFlag_added_left == 0:
                        iFlag_added_left = 1
                        aCoord_gcs_left.append([180-1.0E-8, dLat_bottom])
                        aCoord_gcs_left.append([180-1.0E-8, dLat_top])
                else:
                    pass
            else:
                if i in aIndex:
                    if iFlag_added_left == 0:
                        iFlag_added_left = 1
                        aCoord_gcs_left.append([180-1.0E-8, dLat_bottom])
                        aCoord_gcs_left.append([180-1.0E-8, dLat_top])
                else:
                    aCoord_gcs_left.append(aCoord_gcs[i])

        if dLongitude < 0:
            if i <= np.min(aIndex_dummy):
                aCoord_gcs_right.append(aCoord_gcs[i])
                if i in aIndex:
                    if iFlag_added_right == 0:
                        iFlag_added_right = 1
                        aCoord_gcs_right.append([-180+1.0E-8, dLat_top])
                        aCoord_gcs_right.append([-180+1.0E-8, dLat_bottom])
                    else:
                        pass
            else:
                if i in aIndex:
                    if iFlag_added_right == 0:
                        aCoord_gcs_right.append([-180+1.0E-8, dLat_bottom])
                        aCoord_gcs_right.append([-180+1.0E-8, dLat_top])
                else:
                    aCoord_gcs_right.append(aCoord_gcs[i])

    return [aCoord_gcs_left, aCoord_gcs_right]


def create_cross_idl_polygon(sFilename_geojson):

    pDriver = ogr.GetDriverByName('GeoJSON')
    pDriver = ogr.GetDriverByName('ESRI Shapefile')


    #set up srs
    pSpatialReference = osr.SpatialReference()
    pSpatialReference.ImportFromEPSG(4326)

    if os.path.exists(sFilename_geojson):
        pDriver.DeleteDataSource(sFilename_geojson)

    pDataSource = pDriver.CreateDataSource(sFilename_geojson)
    pLayer = pDataSource.CreateLayer('cross', pSpatialReference, ogr.wkbPolygon)
    #add id field
    pField = ogr.FieldDefn('id', ogr.OFTInteger)
    pLayer.CreateField(pField)

    dLongitude_lu =  -175
    dLongitude_rd =  175
    dLatitude_lu =  5
    dLatitude_rd =  -5

    #if dLongitude_lu < -150 and dLongitude_rd > 150:
    #    dLongitude_lu = dLongitude_lu + 360
    wkt = 'POLYGON ((-179.997342114428 66.2769360301573,179.98848859037 66.2814388694679,179.984323361886 66.281448567623,179.976318838375 66.27728074511330,179.976930616892 66.275689036998,179.98551001155 66.2729876032834,179.990204009696 66.2729766264868,-179.997342114428 66.2769360301573))'

    #create a polygon using wkt
    pGeometry = ogr.CreateGeometryFromWkt(wkt)

    pFeature = ogr.Feature(pLayer.GetLayerDefn())
    #pGeometry = ogr.Geometry(ogr.wkbPolygon)
    #pRing = ogr.Geometry(ogr.wkbLinearRing)
    #pRing.AddPoint(dLongitude_lu, dLatitude_lu)
    #pRing.AddPoint(dLongitude_rd, dLatitude_lu)
    #pRing.AddPoint(dLongitude_rd, dLatitude_rd)
    #pRing.AddPoint(dLongitude_lu, dLatitude_rd)
    #pRing.AddPoint(dLongitude_lu, dLatitude_lu)

    #pGeometry.AddGeometry(pRing)

    if pGeometry.IsValid() == False:
        print('Warning: invalid polygon')

        

    aCoord_gcs = get_geometry_coordinates(pGeometry)

    aCoord_gcs_split = split_polygon_at_idl_manual(aCoord_gcs)

      # Add the split polygons to the layer
    for i in range(2):
        if i == 1:
            continue
        aCoord_gcs = aCoord_gcs_split[i]

        pGeometry = ogr.Geometry(ogr.wkbPolygon)

        pRing = ogr.Geometry(ogr.wkbLinearRing)
        for aCoord in aCoord_gcs:
            pRing.AddPoint(aCoord[0], aCoord[1])

        pRing.AddPoint(aCoord_gcs[0][0], aCoord_gcs[0][1])
        pGeometry.AddGeometry(pRing)

        #get wkt
        pGeometry.FlattenTo2D()
        sWkt = pGeometry.ExportToWkt()
        print(sWkt)
        if pGeometry.IsValid() == False:
            print('Warning: invalid polygon')
        pFeature.SetGeometry(pGeometry)
        #set id
        pFeature.SetField('id', i)
        pLayer.CreateFeature(pFeature)

    #pFeature.SetGeometry(pGeometry)
    #pLayer.CreateFeature(pFeature)

    pDataSource.Destroy()
    #close
    pLayer = None
    pDataSource = None
    pSpatialReference = None




    return

if __name__ == '__main__':

    sFolder_out = '/qfs/people/liao313/workspace/python/hexwatershed_utility/figures/mpas'

    sFilename_geojson = os.path.join(sFolder_out, 'cross_idl.shp')
    create_cross_idl_polygon(sFilename_geojson)