#this function can be used to extract the cell information from the MPAS mesh file and check whether the mesh cell is valid or not
#Author: Chang Liao, chang.liao@pnnl.gov
import os, sys
import math
import importlib.util
import numpy as np
import netCDF4 as nc
from osgeo import gdal, osr, ogr
iFlag_cython = importlib.util.find_spec("cython")
#convert from mpas 360 to 180 degree for longitude
from pyearth.system.define_global_variables import *
iFlag_cython = importlib.util.find_spec("cython")
if iFlag_cython is not None:
    from pyflowline.algorithms.cython.kernel import convert_360_to_180
else:
    from pyearth.gis.geometry.convert_longitude_range import convert_360_to_180
from pyearth.gis.geometry.convert_longitude_range import convert_360_to_180_np
from pyearth.gis.location.xyz_to_lonlat import xyz_to_lonlat
from pyearth.toolbox.data.geoparquet.convert_geojson_to_geoparquet import convert_geojson_to_geoparquet
gdal.UseExceptions()
import math


def exam_jigsaw_mesh(sFilename_jigsaw_mesh_netcdf ,
                             iFlag_plot_in = None           ):

    pDatasets_mesh = nc.Dataset(sFilename_jigsaw_mesh_netcdf, 'r')
    #get netcdf format
    format = pDatasets_mesh.file_format
    print('Netcdf format:', format)
    #read new netcdf
    for sKey, aValue in pDatasets_mesh.variables.items():
        if sKey == 'xCell':
            xCell0 = aValue

        if sKey == 'yCell':
            yCell0 = aValue

        if sKey == 'zCell':
            zCell0 = aValue

        if sKey == 'xVertex':
            xVertex0 = aValue

        if sKey == 'yVertex':
            yVertex0 = aValue

        if sKey == 'zVertex':
            zVertex0 = aValue

        if sKey == 'cellsOnVertex':
            cellOnVertex0 = aValue

    xCell = xCell0[:]
    yCell = yCell0[:]
    zCell = zCell0[:]
    xVertex = xVertex0[:]
    yVertex = yVertex0[:]
    zVertex = zVertex0[:]
    cellOnVertex = cellOnVertex0[:]

    nCell = len(xCell0)
    nVertex = len(xVertex)

    #convert xyz to long  lat
    sFilename_geojson = sFilename_jigsaw_mesh_netcdf.replace('.nc', '_combined.geojson')
    if os.path.exists(sFilename_geojson):
        os.remove(sFilename_geojson)
    pDriver = ogr.GetDriverByName('GeoJSON')
    pDS = pDriver.CreateDataSource(sFilename_geojson)
    #use wgs84
    pSRS_wgs84 = osr.SpatialReference()
    pSRS_wgs84.ImportFromEPSG(4326)
    pLayer = pDS.CreateLayer('point', srs=pSRS_wgs84)
    pFieldDefn = ogr.FieldDefn('id', ogr.OFTInteger)
    pLayer.CreateField(pFieldDefn)
    pFieldDefn = ogr.FieldDefn('type', ogr.OFTInteger)
    pLayer.CreateField(pFieldDefn)
    #add lon and lat to the point
    pFieldDefn = ogr.FieldDefn('lon', ogr.OFTReal)
    pLayer.CreateField(pFieldDefn)
    pFieldDefn = ogr.FieldDefn('lat', ogr.OFTReal)
    pLayer.CreateField(pFieldDefn)

    index = 1
    for i in range(nVertex):
        #pt = Pnt(xVertex[i], yVertex[i], zVertex[i], i)
        #pt.normalize()
        #dLongitude_r =  float(pt.x)
        #dLatitude_r =  float(pt.y)
        #dLongitude360 =  dLongitude_r / math.pi * 180
        #dLongitude = convert_360_to_180(dLongitude360)
        #dLatitude =  dLatitude_r / math.pi * 180

        dLongitude, dLatitude = xyz_to_lonlat(xVertex[i], yVertex[i], zVertex[i])

        pPoint = ogr.Geometry(ogr.wkbPoint)
        pPoint.AddPoint(dLongitude, dLatitude)
        pFeature = ogr.Feature(pLayer.GetLayerDefn())
        pFeature.SetField('id', index)
        pFeature.SetField('type', 1)
        pFeature.SetField('lon', dLongitude)
        pFeature.SetField('lat', dLatitude)
        pFeature.SetGeometry(pPoint)
        pLayer.CreateFeature(pFeature)
        pFeature = None
        index = index + 1

    for i in range(nCell):
        #pt = Pnt(xCell[i], yCell[i], zCell[i], i)
        #pt.normalize()
        #dLongitude_r =  float(pt.x)
        #dLatitude_r =  float(pt.y)
        #dLongitude360 =  dLongitude_r / math.pi * 180
        #dLongitude = convert_360_to_180(dLongitude360)
        #dLatitude =  dLatitude_r / math.pi * 180
        dLongitude, dLatitude = xyz_to_lonlat(xCell[i], yCell[i], zCell[i])

        pPoint = ogr.Geometry(ogr.wkbPoint)
        pPoint.AddPoint(dLongitude, dLatitude)
        pFeature = ogr.Feature(pLayer.GetLayerDefn())
        pFeature.SetField('id', index)
        pFeature.SetField('type', 2)
        pFeature.SetField('lon', dLongitude)
        pFeature.SetField('lat', dLatitude)
        pFeature.SetGeometry(pPoint)
        pLayer.CreateFeature(pFeature)
        pFeature = None
        index = index + 1

    pDS = None
    sFilename_parquet = sFilename_geojson.replace('.geojson', '.parquet')
    convert_geojson_to_geoparquet(sFilename_geojson, sFilename_parquet)

    return
#create a main call
if __name__ == '__main__':

    sFilename_jigsaw_mesh_netcdf = '/qfs/people/liao313/workspace/cplus/mpas_debug/mesh_in.nc'
    iFlag_plot_in = 1
    exam_jigsaw_mesh(sFilename_jigsaw_mesh_netcdf, iFlag_plot_in = 1)
