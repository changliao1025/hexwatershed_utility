#this function injects the DEM into the mesh file, most of the time, the mesh should be the MPAS mesh
import os, sys
import math
import importlib.util
import numpy as np
import netCDF4 as nc
from osgeo import gdal, osr, ogr
from pyearth.toolbox.data.beta.add_variable_to_netcdf import add_variable_to_netcdf

gdal.UseExceptions()

def inject_dem_into_mpas_mesh_netcdf(  sFilename_mesh_netcdf,
                                sFilename_mesh_with_dem_netcdf,
                                sFilename_mesh_geojson ):

    pDriver_geojson = ogr.GetDriverByName('GeoJSON')
    #check file existence
    if not os.path.exists(sFilename_mesh_netcdf):
        print("The mesh file does not exist, please check the file path")
        return

    if not os.path.exists(sFilename_mesh_geojson):
        print("The mesh file does not exist, please check the file path")
        return


    #check output file, remove if exists
    if os.path.exists(sFilename_mesh_with_dem_netcdf):
        os.remove(sFilename_mesh_with_dem_netcdf)


    #read mesh geojson file
    pDataSource = pDriver_geojson.Open(sFilename_mesh_geojson, 0)
    pLayer = pDataSource.GetLayer()
    iFeatureCount = pLayer.GetFeatureCount()
    #only need the field value of all the features
    aBed_elevation = np.zeros(iFeatureCount)
    aIce_thickness = np.zeros(iFeatureCount)
    pFeature = pLayer.GetNextFeature()
    iIndex = 0
    while pFeature:
        fBed_elevation = pFeature.GetField('bed_elevation')
        fIce_thickness = pFeature.GetField('ice_thickness')
        aBed_elevation[iIndex] = float(fBed_elevation)
        aIce_thickness[iIndex] = float(fIce_thickness)
        pFeature = pLayer.GetNextFeature()
        iIndex = iIndex + 1

    #convert to numpy array
    aBed_elevation = np.array(aBed_elevation)
    aIce_thickness = np.array(aIce_thickness)
    #add into the netcdf file

    sFilename_old=sFilename_mesh_netcdf
    sFilename_new=sFilename_mesh_with_dem_netcdf
    aData_in= [aBed_elevation, aIce_thickness]
    sVariable_in= ['bed_elevation', 'ice_thickness']
    sUnit_in= ['m', 'm']
    aDimension_in = [['nCells'], ['nCells']]
    add_variable_to_netcdf(sFilename_old, sFilename_new, aData_in, sVariable_in, sUnit_in, aDimension_in)

    return

if __name__ == '__main__':
    sFilename_mesh_netcdf = '/compyfs/liao313/04model/pyflowline/conus/pyflowline20241201019/jigsaw/out/invert_mesh.nc'
    sFilename_mesh_with_dem_netcdf='/compyfs/liao313/00raw/mesh/global/hrm/mpas_with_land_elevation.nc'
    sFilename_mesh_geojson = '/compyfs/liao313/00raw/mesh/global/hrm/mpas_land_with_elevation.geojson'

    inject_dem_into_mpas_mesh_netcdf(  sFilename_mesh_netcdf,
                                sFilename_mesh_with_dem_netcdf,
                                sFilename_mesh_geojson    )