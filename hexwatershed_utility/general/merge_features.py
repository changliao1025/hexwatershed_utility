import os, sys
import numpy as np
from osgeo import gdal, osr, ogr, gdalconst

from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates
from pyearth.gis.geometry.calculate_polygon_area import calculate_polygon_area
def merge_features(sFilename_in, sFilename_out, sFormat='GeoJSON', iFlag_exterior_in = 0):

    pDataset_in = ogr.Open(sFilename_in)
    # Get the first layer in the file
    pLayer_in = pDataset_in.GetLayer(0)
    # Count the number of features (polygons)
    nFeature = pLayer_in.GetFeatureCount()
    # Get the spatial reference of the layer
    pSpatial_reference = pLayer_in.GetSpatialRef()
    wkt2 = pSpatial_reference.ExportToWkt()

    if nFeature == 0:
        print('No feature found in the input shapefile')
        return
    else:
        print('Number of features found in the input shapefile: ' + str(nFeature))

    # Create a new dataset using the output filename
    pDriver = ogr.GetDriverByName(sFormat)
    if pDriver is None:
        print('Driver not found')
        return

    if os.path.exists(sFilename_out):
        pDriver.DeleteDataSource(sFilename_out)

    #obtain the geotype of first layer and
    iGeomType = pLayer_in.GetGeomType()
    #obtain the geotype of first geometry
    pLayer_in.ResetReading()

    # Obtain the first feature
    pFeature = pLayer_in.GetNextFeature()
    pGeometry = pFeature.GetGeometryRef()
    #flatten to 2D
    pGeometry.FlattenTo2D()
    if pGeometry is None:
        print('Geometry not found')
        return
    #get the geometry type
    iGeomType = pGeometry.GetGeometryType()
    if iGeomType == ogr.wkbPoint or iGeomType == ogr.wkbLineString or iGeomType == ogr.wkbPolygon:
        #check whether it is a multi-geometry
        pass
    else:
        if iGeomType == ogr.wkbMultiPoint or iGeomType == ogr.wkbMultiLineString or iGeomType == ogr.wkbMultiPolygon:
            #get the number of geometries
            nGeom = pGeometry.GetGeometryCount()
            #get the first geometry
            pGeometry_single = pGeometry.GetGeometryRef(0)
            iGeomType = pGeometry_single.GetGeometryType()
        else:
            sGeomType = ogr.GeometryTypeToName(iGeomType)
            print('Geometry type not supported', sGeomType)

    #create a temp file if exterior is needed
    if iFlag_exterior_in == 1:
        sFilename_out_temp = sFilename_out + '_temp'
    else:
        sFilename_out_temp = sFilename_out

    if os.path.exists(sFilename_out_temp):
        pDriver.DeleteDataSource(sFilename_out_temp)

    pDataset_out = pDriver.CreateDataSource(sFilename_out_temp)
    if pDataset_out is None:
        print('Dataset not created')
        return

    if iGeomType == ogr.wkbPoint:
        #create the layer
        pLayer_out = pDataset_out.CreateLayer('layer', pSpatial_reference, geom_type=ogr.wkbPoint)
        pGeometry_merge = ogr.Geometry(ogr.wkbPoint)
    else:
        if iGeomType == ogr.wkbLineString:
            #create the layer
            pLayer_out = pDataset_out.CreateLayer('layer', pSpatial_reference, geom_type=ogr.wkbLineString)
            pGeometry_merge = ogr.Geometry(ogr.wkbLineString)
        else:
            if iGeomType == ogr.wkbPolygon:
                #create the layer
                pLayer_out = pDataset_out.CreateLayer('layer', pSpatial_reference, geom_type=ogr.wkbPolygon)
                pGeometry_merge = ogr.Geometry(ogr.wkbPolygon)
            else:
                #convert to name
                sGeomType = ogr.GeometryTypeToName(iGeomType)
                print('Geometry type not supported', sGeomType)
                return
            pass
        pass
    # Create a new layer in the output shapefile

    # Loop through the input features and merge them into the output layer
    pLayer_in.ResetReading()  # Reset reading to start from the first feature again
    pFeature = pLayer_in.GetNextFeature()
    while pFeature:
        pGeometry = pFeature.GetGeometryRef()
        pGeometry.FlattenTo2D()
        if pGeometry is not None:
            #check geotype again
            iGeomType_new = pGeometry.GetGeometryType()
            if iGeomType_new == iGeomType:
                # Union the geometry of each feature with the merged polygon
                pGeometry_merge = pGeometry_merge.Union(pGeometry)
            else:
                #check whether the geometry type is a multi-geometry
                if iGeomType_new == ogr.wkbMultiPoint or iGeomType_new == ogr.wkbMultiLineString or iGeomType_new == ogr.wkbMultiPolygon:
                    #get the number of geometries
                    nGeom = pGeometry.GetGeometryCount()
                    aArea = np.zeros(nGeom)
                    for i in range(nGeom):
                        pGeometry_single = pGeometry.GetGeometryRef(i)
                        #check again its geometry type
                        iGeomType_single = pGeometry_single.GetGeometryType()
                        if iGeomType_single == iGeomType:
                            #calculate the area
                            aCoords_gcs = get_geometry_coordinates(pGeometry_single)
                            dArea = calculate_polygon_area(aCoords_gcs[:,0], aCoords_gcs[:,1])
                            aArea[i] = dArea

                    #get the index of the largest area
                    iIndex = np.argmax(aArea)
                    #get the geometry
                    pGeometry_single = pGeometry.GetGeometryRef(int(iIndex))
                    # Union the geometry of each feature with the merged polygon
                    pGeometry_merge = pGeometry_merge.Union(pGeometry_single)
                else:
                    #convert to name
                    sGeomType = ogr.GeometryTypeToName(iGeomType_new)
                    print('Geometry type not supported')

            pFeature = pLayer_in.GetNextFeature()

    # Create a new feature in the output layer
    pFeature_out = ogr.Feature(pLayer_out.GetLayerDefn())
    pFeature_out.SetGeometry(pGeometry_merge)
    pLayer_out.CreateFeature(pFeature_out)

    #close the dataset
    pDataset_out.Destroy()
    pDataset_in.Destroy()

    if iFlag_exterior_in == 1:
        #now we will use the temp file to obtain the correct geometry
        pDataset_out = ogr.Open(sFilename_out_temp)
        pLayer_out = pDataset_out.GetLayer(0)
        pFeature_out = pLayer_out.GetNextFeature()
        pGeometry_out = pFeature_out.GetGeometryRef()
        #check geotype again
        iGeomType_new = pGeometry_out.GetGeometryType()
        sGeomType = ogr.GeometryTypeToName(iGeomType_new)
        print('Geometry type of the merged feature: ' , sGeomType)
        #get the external ring of the polygon
        aCoords_exterior = []
        pExteriorRing = pGeometry_out.GetGeometryRef(0)
        #get the number of points
        nPoints = pExteriorRing.GetPointCount()
        for i in range(nPoints):
            point = pExteriorRing.GetPoint(i)
            aCoords_exterior.append((point[0], point[1]))
        #save the external ring as a polygon
        pGeometry_out = ogr.Geometry(ogr.wkbPolygon)
        ring = ogr.Geometry(ogr.wkbLinearRing)
        for coord in aCoords_exterior:
            ring.AddPoint(coord[0], coord[1])
        pGeometry_out.AddGeometry(ring)

        #use SetPricision to remove the Slivers
        pGeometry_out.SetPrecision(8)
        pDataset_out = pDriver.CreateDataSource(sFilename_out)
        if pDataset_out is None:
            print('Dataset not created')
            return
        pLayer_out = pDataset_out.CreateLayer('layer', pSpatial_reference, geom_type=ogr.wkbPolygon)
        pFeature_out = ogr.Feature(pLayer_out.GetLayerDefn())
        #dBuffer_threshold_in = -0.0001
        #pGeometry_out = pGeometry_out.Buffer(dBuffer_threshold_in)
        #pGeometry_out = pGeometry_out.Buffer(-dBuffer_threshold_in)
        pFeature_out.SetGeometry(pGeometry_out)
        pLayer_out.CreateFeature(pFeature_out)
        pDataset_out.Destroy()

        print('Finished merging features.')

        #delete the temp file
        pDriver.DeleteDataSource(sFilename_out_temp)

    return