#this function injects the DEM into the mesh file, most of the time, the mesh should be the MPAS mesh
import os, sys
import math
import importlib.util
import contextlib
import numpy as np
import netCDF4 as nc
from osgeo import gdal, osr, ogr
gdal.UseExceptions()
gdal.PushErrorHandler('CPLQuietErrorHandler')
iFlag_cython = importlib.util.find_spec("cython")
if iFlag_cython is not None:
    from tinyr import RTree

from pyearth.gis.gdal.gdal_vector_format_support import get_vector_driver_from_filename
from pyearth.gis.gdal.gdal_validate_polygon_file import gdal_validate_polygon_file
from pyearth.gis.geometry.convert_longitude_range import convert_360_to_180_np

def convert_mpas_mesh_to_vector( sFilename_mpas_mesh_netcdf_base,
                                 sFilename_vector_out ):

    pDriver_out = get_vector_driver_from_filename(sFilename_vector_out)
    #check file existence
    if not os.path.exists(sFilename_mpas_mesh_netcdf_base):
        print("The mesh file does not exist, please check the file path")
        return

    #check output file, remove if exists

    if os.path.exists(sFilename_vector_out):
        os.remove(sFilename_vector_out)



    #first, we need to read the mesh file and copy it to the new file
    pDatasets_mesh = nc.Dataset(sFilename_mpas_mesh_netcdf_base, 'r')
    #get netcdf format
    format = pDatasets_mesh.file_format

    #read new netcdf
    for sKey, aValue in pDatasets_mesh.variables.items():
        #we need to filter out unused grids based on mpas specs
        if sKey == 'lonCell':
            lonCell0 = aValue

        if sKey == 'latCell':
            latCell0 = aValue

        if sKey == 'verticesOnCell':
            verticesOnCell0 = aValue

        if sKey == 'indexToCellID':
            indexToCellID0 = aValue

        if sKey == 'lonVertex':
            lonVertex0 = aValue

        if sKey == 'latVertex':
            latVertex0 = aValue

    aLongitudeCell = lonCell0[:] / math.pi * 180
    aLatitudeCell = latCell0[:] / math.pi * 180
    aLatitudeVertex = latVertex0[:] / math.pi * 180
    aLongitudeVertex = lonVertex0[:] / math.pi * 180
    aVertexOnCell = verticesOnCell0[:]
    aIndexToCellID = indexToCellID0[:]
    ncell = len(aIndexToCellID)

    pSpatialRef_target = osr.SpatialReference()
    pSpatialRef_target.ImportFromEPSG(4326)

    pDataset = pDriver_out.CreateDataSource(sFilename_vector_out)
    pLayerOut = pDataset.CreateLayer('cell', pSpatialRef_target, ogr.wkbPolygon)
    #create field for id, lon, lat
    pFieldDefn = ogr.FieldDefn('id', ogr.OFTInteger)
    pLayerOut.CreateField(pFieldDefn)
    pFieldDefn = ogr.FieldDefn('lon', ogr.OFTReal)
    #set width for lon and lat
    pFieldDefn.SetWidth(25)
    pFieldDefn.SetPrecision(15)
    pLayerOut.CreateField(pFieldDefn)
    pFieldDefn = ogr.FieldDefn('lat', ogr.OFTReal)
    pFieldDefn.SetWidth(25)
    pFieldDefn.SetPrecision(15)
    pLayerOut.CreateField(pFieldDefn)
    pLayerDefn = pLayerOut.GetLayerDefn()
    aLongitudeCell_180 = convert_360_to_180_np(aLongitudeCell)
    aLongitudeVertex_180 = convert_360_to_180_np(aLongitudeVertex)

    for i in range(ncell):
        dLon_center = float(aLongitudeCell_180[i])
        dLat_center = float(aLatitudeCell[i])
        lCellID = int(aIndexToCellID[i])
        aVertexOnCellIndex = np.array(aVertexOnCell[i,:])
        dummy0 = np.where(aVertexOnCellIndex > 0)
        aVertexIndex = aVertexOnCellIndex[dummy0] - 1
        aLonVertex = aLongitudeVertex_180[aVertexIndex]
        aLatVertex = aLatitudeVertex[aVertexIndex]
        nVertex = len(aLonVertex)
        if nVertex < 3:
            print("Vertex number is: ", nVertex, i)
        else:
            ring = ogr.Geometry(ogr.wkbLinearRing)
            for lon, lat in zip(aLonVertex, aLatVertex):
                ring.AddPoint(float(lon), float(lat))

            ring.CloseRings()
            pPolygon = ogr.Geometry(ogr.wkbPolygon)
            pPolygon.AddGeometry(ring)
            if pPolygon.IsValid():
                pFeatureOut = ogr.Feature(pLayerDefn)
                pFeatureOut.SetGeometry(pPolygon)
                pFeatureOut.SetField('id', lCellID)
                pFeatureOut.SetField('lon', dLon_center)
                pFeatureOut.SetField('lat', dLat_center)
                pLayerOut.CreateFeature(pFeatureOut)

    #close the dataset
    pDataset.FlushCache()
    pDatasets_mesh.close()
    pFeatureOut = None
    pLayerOut = None
    pDataset = None

    #check whether all the polygons are valid
    iFlag_valid = gdal_validate_polygon_file(sFilename_vector_out)
    print("The validity of the polygon is: ", iFlag_valid)

    return
