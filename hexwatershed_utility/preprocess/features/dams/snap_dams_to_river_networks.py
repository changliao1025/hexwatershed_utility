import os, sys
import numpy as np
from osgeo import ogr, osr
from rtree.index import Index as RTreeindex
from pyflowline.classes.vertex import pyvertex
from pyflowline.formats.read_flowline import read_flowline_geojson
from pyearth.gis.gdal.gdal_vector_format_support import get_vector_driver_from_filename

def snap_dams_to_river_networks(sFilename_dam, sFilename_river, sFilename_out, sVariable_drainage='drainage'):
    if not os.path.exists(sFilename_dam):
        print('dam location file does not exist')
        sys.exit()

    if not os.path.exists(sFilename_river):
        print('river network does not exist')
        sys.exit()


    pDataSource_dam = ogr.Open(sFilename_dam, 0)
    pLayer_dam = pDataSource_dam.GetLayer()
    nDam = pLayer_dam.GetFeatureCount()

    pSpatialReference = osr.SpatialReference()
    pSpatialReference.ImportFromEPSG(4326)

    # Create three output files for different scenarios
    # Extract file extension and base name
    sFilename_base, sExt = os.path.splitext(sFilename_out)
    sFilename_on = f"{sFilename_base}_on{sExt}"
    sFilename_snapped = f"{sFilename_base}_snapped{sExt}"
    sFilename_merged = f"{sFilename_base}_merged{sExt}"
    sFilename_no_river = f"{sFilename_base}_no_river{sExt}"
    sFilename_too_far = f"{sFilename_base}_too_far{sExt}"

    pDriver = get_vector_driver_from_filename(sFilename_out)

    # Remove existing files if they exist
    for sFile in [sFilename_on, sFilename_snapped, sFilename_no_river, sFilename_too_far]:
        if os.path.exists(sFile):
            pDriver.DeleteDataSource(sFile)

    # Create three data sources and layers
    pDataSource_on = pDriver.CreateDataSource(sFilename_on)
    pLayer_on = pDataSource_on.CreateLayer('dam', pSpatialReference, ogr.wkbPoint)

    pDataSource_snapped = pDriver.CreateDataSource(sFilename_snapped)
    pLayer_snapped = pDataSource_snapped.CreateLayer('dam', pSpatialReference, ogr.wkbPoint)

    pDataSource_no_river = pDriver.CreateDataSource(sFilename_no_river)
    pLayer_no_river = pDataSource_no_river.CreateLayer('dam', pSpatialReference, ogr.wkbPoint)

    pDataSource_too_far = pDriver.CreateDataSource(sFilename_too_far)
    pLayer_too_far = pDataSource_too_far.CreateLayer('dam', pSpatialReference, ogr.wkbPoint)

    pDataSource_merged = pDriver.CreateDataSource(sFilename_merged)
    pLayer_merged = pDataSource_merged.CreateLayer('dam', pSpatialReference, ogr.wkbPoint)

    # Create fields for all three layers
    for pLayer in [pLayer_snapped, pLayer_no_river, pLayer_too_far, pLayer_on, pLayer_merged]:
        pLayer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
        pLayer.CreateField(ogr.FieldDefn('dLongitude_degree', ogr.OFTReal))
        pLayer.CreateField(ogr.FieldDefn('dLatitude_degree', ogr.OFTReal))
        pLayer.CreateField(ogr.FieldDefn('drainage', ogr.OFTReal))

    pLayerDefn_snapped = pLayer_snapped.GetLayerDefn()
    pLayerDefn_no_river = pLayer_no_river.GetLayerDefn()
    pLayerDefn_too_far = pLayer_too_far.GetLayerDefn()
    pLayerDefn_on = pLayer_on.GetLayerDefn()
    pLayerDefn_merged = pLayer_merged.GetLayerDefn()

    pDataSource_river = ogr.Open(sFilename_river, 0)
    pLayer_river = pDataSource_river.GetLayer()
    nFlowline = pLayer_river.GetFeatureCount()

    aFlowline, pProjection_geojson = read_flowline_geojson( sFilename_river )

    index_flowline = RTreeindex()
    for i in range(nFlowline):
        pBound= aFlowline[i].pBound
        index_flowline.insert(i, pBound)  #
        pass

    #read dam location
    iDam_on = 1
    iDam_snapped = 1
    iDam_no_river = 1
    iDam_too_far = 1
    nCount_on = 0
    nCount_snapped = 0
    nCount_no_river = 0
    nCount_too_far = 0
    dBuffer = 0.01 #degree, approximately 1 km

    for i in range(nDam):
        pFeature_dam = pLayer_dam.GetFeature(i)
        point = pFeature_dam.GetGeometryRef()
        dDrainage = pFeature_dam.GetField(sVariable_drainage)
        #get point coordinates using gdal api
        x = point.GetX()
        y = point.GetY()
        point0= dict()
        point0['dLongitude_degree'] = x
        point0['dLatitude_degree'] = y
        pVertex=pyvertex(point0)

        left =   x - dBuffer
        right =  x + dBuffer
        bottom = y - dBuffer
        top =    y + dBuffer
        pBound= (left, bottom, right, top)
        aIntersect = list(index_flowline.intersection( pBound )  )

        if len(aIntersect) == 0:
            # Scenario 2: Dam is not close to any river
            pPoint = ogr.Geometry(ogr.wkbPoint)
            pPoint.AddPoint(x, y)
            pFeature_out = ogr.Feature(pLayerDefn_no_river)
            pFeature_out.SetGeometry(pPoint)
            pFeature_out.SetField('id', iDam_no_river)
            pFeature_out.SetField('dLongitude_degree', x)
            pFeature_out.SetField('dLatitude_degree', y)
            pFeature_out.SetField('drainage', dDrainage)
            pLayer_no_river.CreateFeature(pFeature_out)
            iDam_no_river = iDam_no_river + 1
            nCount_no_river = nCount_no_river + 1
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
                if dDistance_min< 10.0:
                    #the dam is already on the river network
                    pPoint = ogr.Geometry(ogr.wkbPoint)
                    pPoint.AddPoint(pVertex_closest.dLongitude_degree, pVertex_closest.dLatitude_degree)
                    pFeature_out = ogr.Feature(pLayerDefn_on)
                    pFeature_out.SetGeometry(pPoint)
                    pFeature_out.SetField('id', iDam_on)
                    pFeature_out.SetField('dLongitude_degree', pVertex_closest.dLongitude_degree)
                    pFeature_out.SetField('dLatitude_degree', pVertex_closest.dLatitude_degree)
                    pFeature_out.SetField('drainage', dDrainage)
                    pLayer_on.CreateFeature(pFeature_out)
                    iDam_on = iDam_on + 1
                    nCount_on = nCount_on + 1
                else:
                # Scenario 1: Dam successfully snapped to river network
                    pPoint = ogr.Geometry(ogr.wkbPoint)
                    pPoint.AddPoint(pVertex_closest.dLongitude_degree, pVertex_closest.dLatitude_degree)
                    pFeature_out = ogr.Feature(pLayerDefn_snapped)
                    pFeature_out.SetGeometry(pPoint)
                    pFeature_out.SetField('id', iDam_snapped)
                    pFeature_out.SetField('dLongitude_degree', pVertex_closest.dLongitude_degree)
                    pFeature_out.SetField('dLatitude_degree', pVertex_closest.dLatitude_degree)
                    pFeature_out.SetField('drainage', dDrainage)
                    pLayer_snapped.CreateFeature(pFeature_out)
                    iDam_snapped = iDam_snapped + 1
                    nCount_snapped = nCount_snapped + 1

                #merge as well
                pPoint = ogr.Geometry(ogr.wkbPoint)
                pPoint.AddPoint(pVertex_closest.dLongitude_degree, pVertex_closest.dLatitude_degree)
                pFeature_out = ogr.Feature(pLayerDefn_merged)
                pFeature_out.SetGeometry(pPoint)
                pFeature_out.SetField('id', i + 1)
                pFeature_out.SetField('dLongitude_degree', pVertex_closest.dLongitude_degree)
                pFeature_out.SetField('dLatitude_degree', pVertex_closest.dLatitude_degree)
                pFeature_out.SetField('drainage', dDrainage)
                pLayer_merged.CreateFeature(pFeature_out)

            else:
                # Scenario 3: Dam is too far from any river (distance >= 3km)
                print('dam is too far from any river')
                print(i, x, y, 'distance is ', dDistance_min)
                pPoint = ogr.Geometry(ogr.wkbPoint)
                pPoint.AddPoint(x, y)
                pFeature_out = ogr.Feature(pLayerDefn_too_far)
                pFeature_out.SetGeometry(pPoint)
                pFeature_out.SetField('id', iDam_too_far)
                pFeature_out.SetField('dLongitude_degree', x)
                pFeature_out.SetField('dLatitude_degree', y)
                pFeature_out.SetField('drainage', dDrainage)
                pLayer_too_far.CreateFeature(pFeature_out)
                iDam_too_far = iDam_too_far + 1
                nCount_too_far = nCount_too_far + 1

    pDataSource_dam.Destroy()
    pDataSource_river.Destroy()
    pDataSource_snapped.Destroy()
    pDataSource_no_river.Destroy()
    pDataSource_on.Destroy()
    pDataSource_too_far.Destroy()

    print('Finished!')
    print(f'Summary:')
    print(f'  Total dams processed: {nDam}')
    print(f'  Already on river network (<10m): {nCount_on} -> {sFilename_on}')
    print(f'  Successfully snapped (< 3km): {nCount_snapped} -> {sFilename_snapped}')
    print(f'  No river nearby: {nCount_no_river} -> {sFilename_no_river}')
    print(f'  Too far from river (>= 3km): {nCount_too_far} -> {sFilename_too_far}')
