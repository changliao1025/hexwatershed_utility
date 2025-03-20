import os, sys
import copy
import numpy as np
import importlib.util
import netCDF4 as nc
from osgeo import osr, ogr, gdal
from pyearth.system.define_global_variables import *
iFlag_cython = importlib.util.find_spec("cython")

if iFlag_cython is not None:
    from tinyr import RTree
    iFlag_use_rtree = 1

def find_gage_mesh_cell_id(aSitename_in, aLongitude_gage_in, aLatitude_gage_in, aDrainage_area_in,
                            sFilename_domain_in,
                            sFilename_parameter_in,
                              iSearch_radius_in = None,
                              dThreshold_drainage_in = None,
                              dThreshold_difference_in = 0.10,
                              iFlag_data_km_in = 0,
                              dBuffer_in = 0.3):

    if dThreshold_drainage_in is not None:
        dThreshold_drainage = dThreshold_drainage_in / 1.0E6
    else:
        dThreshold_drainage = 1.0E2
    # Find the cell ID of the gage location using spatial index search
    sFilename_geojson_mesh = sFilename_parameter_in.replace('.nc', '_mesh.geojson')
    if os.path.exists(sFilename_geojson_mesh):
        os.remove(sFilename_geojson_mesh)

    pDatasets_domain = nc.Dataset(sFilename_domain_in, 'r')
    pDimension = pDatasets_domain.dimensions.keys()
    for sKey, aValue in pDatasets_domain.variables.items():
        if (sKey == 'xv'):
            aXV = (aValue[:]).data
            continue
        if (sKey == 'yv'):
            aYV = (aValue[:]).data
            continue
        if (sKey == 'xc'):
            aXC = (aValue[:]).data
            continue
    print(sFilename_parameter_in)
    aDatasets = nc.Dataset(sFilename_parameter_in)
    netcdf_format = aDatasets.file_format
    # Copy variables

    if iFlag_use_rtree == 1:
        #read mesh using the
        index_cell = RTree(max_cap=5, min_cap=2)
    else:
        print('Rtree is not supported')

    for sKey, aValue in aDatasets.variables.items():
        #print(sKey, aValue)
        print(aValue.datatype)
        print(aValue.dimensions)
        # Copy variable attributes
        #outVar.setncatts({k: aValue.getncattr(k) for k in aValue.ncattrs()})
        if sKey == 'CellID':
            aCellID =  (aValue[:]).data
            iFlag_global_id = 1
        if sKey == 'ID':
            aID =  (aValue[:]).data
        if sKey == 'dnID':
            aDnID =  (aValue[:]).data
        if sKey == 'fdir':
            aFdir =  (aValue[:]).data
        if sKey == 'latixy':
            aLatitude = (aValue[:]).data
        if sKey == 'longxy':
            aLongitude = (aValue[:]).data
        if sKey == 'areaTotal2':
            if iFlag_data_km_in == 1:
                aAccu = (aValue[:]).data
            else:
                aAccu = (aValue[:]).data  / 1.0e+6 #from m2 to km2

    pDriver = ogr.GetDriverByName('GeoJSON')
    pDataset = pDriver.CreateDataSource(sFilename_geojson_mesh)
    pSpatial_reference_gcs = osr.SpatialReference()
    pSpatial_reference_gcs.ImportFromEPSG(4326)    # WGS84 lat/long
    #pLayer = pDataset.CreateLayer(sVariable_parameter, pSrs, ogr.wkbPoint)
    pLayer = pDataset.CreateLayer('cell', pSpatial_reference_gcs, ogr.wkbPolygon)
    pLayer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger64))
    pLayer.CreateField(ogr.FieldDefn('cellid', ogr.OFTInteger64))
    pLayer.CreateField(ogr.FieldDefn('dnID', ogr.OFTInteger64))
    pLayer.CreateField(ogr.FieldDefn('drain', ogr.OFTReal))

    pLayerDefn = pLayer.GetLayerDefn()
    pFeature = ogr.Feature(pLayerDefn)
    nCell = len(aID)
    aBoundary = list()
    #also create a large boundary that combines all the boundaries
    #call gdal to create a large boundary geometry
    pGeometry_union  = ogr.Geometry(ogr.wkbPolygon)
    for i in np.arange(nCell):
        ring = ogr.Geometry(ogr.wkbLinearRing)
        pXV0 = aXV[i]
        #remove the dummy value
        pXV = pXV0[pXV0 != -9999]
        nVertex = len(pXV)
        for j in range(nVertex):
            x1 = aXV[i,0,j ]
            y1 = aYV[i,0,j]
            ring.AddPoint(x1, y1)
            pass
        #add the first point to close the polygon
        ring.AddPoint(aXV[i,0,0], aYV[i,0,0])
        pPolygon = ogr.Geometry(ogr.wkbPolygon)
        pPolygon.AddGeometry(ring)
        aBoundary.append(pPolygon)
        pFeature.SetGeometry(pPolygon)
        lID = aID[i]
        #lCellID= aCellID[i]
        lID_down = aDnID[i]
        #if(lID_down != -9999):
            #define id first
        pFeature.SetField('id', lID )
        #pFeature.SetField('cellid', lCellID )
        pFeature.SetField( 'dnID', aDnID[i] )
        pFeature.SetField( 'drain', aAccu[i] )
        #now create the feature
        pLayer.CreateFeature(pFeature)
        #add into tree
        left,  right, bottom,top = pPolygon.GetEnvelope()
        pBound= (left, bottom, right, top)
        index_cell.insert(i, pBound)  #

        pGeometry_union = pGeometry_union.Union(pPolygon)

    #Save and close everything
    pDataset = pLayer = pFeature  = None

    #then process each gaga
    nGage  = len(aLongitude_gage_in)

    #also create a point based vector
    sFilename_geojson_gage = sFilename_parameter_in.replace('.nc', '_gage.geojson')
    if os.path.exists(sFilename_geojson_gage):
        os.remove(sFilename_geojson_gage)
    pDataset_gage = pDriver.CreateDataSource(sFilename_geojson_gage)
    pLayer_gage = pDataset_gage.CreateLayer('gage', pSpatial_reference_gcs, ogr.wkbPoint)
    pLayer_gage.CreateField(ogr.FieldDefn('name', ogr.OFTString))
    pLayer_gage.CreateField(ogr.FieldDefn('lon', ogr.OFTReal))
    pLayer_gage.CreateField(ogr.FieldDefn('lat', ogr.OFTReal))
    pLayer_gage.CreateField(ogr.FieldDefn('drain', ogr.OFTReal))
    pLayer_gage.CreateField(ogr.FieldDefn('drain1', ogr.OFTReal))
    pLayer_gage.CreateField(ogr.FieldDefn('qa', ogr.OFTReal))
    pLayer_gage.CreateField(ogr.FieldDefn('cellid', ogr.OFTInteger64))
    pLayerDefn = pLayer_gage.GetLayerDefn()

    aIndex_out  = list()
    aCellID_out = list()
    aDrainage_area_out = list()

    for i in range(nGage):
        sName = aSitename_in[i]
        dLongitude_gage = float(aLongitude_gage_in[i])
        dLatitude_gage = float(aLatitude_gage_in[i])
        dDrainage_area = float(aDrainage_area_in[i])

        #there is a chance the gsim gage has no drainage area
        #check whether dDrainage_area is nan
        if np.isnan(dDrainage_area):
            iFlag_drainage = 0
            continue
        #create a point
        pPoint = ogr.Geometry(ogr.wkbPoint)
        pPoint.AddPoint(dLongitude_gage, dLatitude_gage)
        if pPoint.Within(pGeometry_union): #check if it is inside the boundary
            pass
        else:
            continue
        #this function only seach a point within a polygon, but it may miss some points
        #aIntersect = list(index_cell.search_surrounding([dLongitude_gage, dLatitude_gage]))
        #use the new search method to find the cell, it uses a rectangle instead of point
        #we will search in a 6*6 matrix, the resolution is 10km, so the buffer is 10 * 3=30km, convert to degree is close to 0.3

        dBuffer = dBuffer_in
        left =   dLongitude_gage - dBuffer
        right =  dLongitude_gage + dBuffer
        bottom = dLatitude_gage - dBuffer
        top =    dLatitude_gage + dBuffer
        pBound= (left, bottom, right, top)

        aIntersect = list(index_cell.search( pBound )  )

        iFlag_qa = 0
        if len(aIntersect) == 0:
            #print('No cell found for the gage')
            continue
        else:
            nIntersect = len(aIntersect)
            #print(nIntersect)
            #set diff as infinity large
            dDiff_min = np.inf
            for j in range(nIntersect):
                lIndex = int(aIntersect[j])
                dDrainage_mesh_cell = float(aAccu[lIndex])
                lDownID = int(aDnID[lIndex])
                dDiff_new = float(np.abs(dDrainage_mesh_cell - dDrainage_area))

                if lDownID < 0:
                    continue

                if dDiff_new < dDiff_min:
                    dDiff_min = float(dDiff_new)
                    dIndex_min = lIndex
                    dDrainage_mesh_cell_min = dDrainage_mesh_cell
            #

            #iFlag = False
            pBoundary=aBoundary[dIndex_min]
            lCellID = aID[dIndex_min]
            #if pPoint.Within(pBoundary):
            #    iFlag = True
                #compare the drainage area difference
                #print(i, dDrainage_mesh_cell, dDrainage_area)
            dRatio = float(dDiff_min / np.max([dDrainage_area, dDrainage_mesh_cell_min]))
            iFlag_qa = float(dRatio)
            if dRatio < dThreshold_difference_in and dDrainage_area > dThreshold_drainage:
                #create a point
                aIndex_out.append(i)
                aCellID_out.append(lCellID)
                aDrainage_area_out.append(dDrainage_mesh_cell_min)
            else:
                print('Warning: the drainage area difference is larger than ', dThreshold_difference_in)

            pFeature = ogr.Feature(pLayerDefn)
            pFeature.SetGeometry(pPoint)
            pFeature.SetField('name', sName)
            pFeature.SetField('lon', dLongitude_gage)
            pFeature.SetField('lat', dLatitude_gage)
            pFeature.SetField('drain', dDrainage_area)
            pFeature.SetField('drain1', dDrainage_mesh_cell_min)
            pFeature.SetField('qa', iFlag_qa)
            pFeature.SetField('cellid', lCellID)
            pLayer_gage.CreateFeature(pFeature)


    pDataset_gage = pLayer_gage = pFeature = None
    return aIndex_out, aCellID_out, aDrainage_area_out