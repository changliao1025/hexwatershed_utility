import os, sys
from osgeo import gdal, ogr, osr
import numpy as np
from pyearth.system.define_global_variables import *

def gdal_warp_antimeridian_polygon(sFilename_shapefile_cut, aFilename_geotiff_in):

    #use the existing polygon
    pDriver_shapefile = ogr.GetDriverByName('ESRI Shapefile')

    pDataset_shapefile = pDriver_shapefile.Open(sFilename_shapefile_cut)
    pLayer_shapefile = pDataset_shapefile.GetLayer()
    pFeature_shapefile = pLayer_shapefile.GetNextFeature()
    pPolygon = pFeature_shapefile.GetGeometryRef()

    pDataset = gdal.Open(aFilename_geotiff_in[0])
    pProjection = pDataset.GetProjection()
    pGeoTransform = pDataset.GetGeoTransform()
    dPixelWidth = pGeoTransform[1]
    pPixelHeight = pGeoTransform[5]
    pBand = pDataset.GetRasterBand(1)
    dMissing_value = pBand.GetNoDataValue()

    pSpatialRef_target = osr.SpatialReference()
    pSpatialRef_target.ImportFromWkt(pProjection)
    pDriver = gdal.GetDriverByName('MEM')


    sResampleAlg = 'near'


    pDataset3 = pDriver_shapefile.CreateDataSource(sFilename_shapefile_cut)
    pLayerOut3 = pDataset3.CreateLayer('cell', pSpatialRef_target, ogr.wkbPolygon)
    pLayerDefn3 = pLayerOut3.GetLayerDefn()
    pFeatureOut3 = ogr.Feature(pLayerDefn3)
    pFeatureOut3.SetGeometry(pPolygon)
    pLayerOut3.CreateFeature(pFeatureOut3)
    pDataset3.FlushCache()
    pWrapOption = gdal.WarpOptions( cropToCutline=True,
                                   cutlineDSName = sFilename_shapefile_cut ,#could be true if vector file is provided
                            xRes=dPixelWidth,
                           yRes=abs(pPixelHeight),
                                dstSRS=pSpatialRef_target , format = 'MEM',
                                resampleAlg=sResampleAlg )
    pDataset_clip_warped = gdal.Warp('', aFilename_geotiff_in, options=pWrapOption)
    aData_clip = pDataset_clip_warped.ReadAsArray()
    aData_clip[aData_clip == dMissing_value] = -9999
    dummy_index0 = np.where(aData_clip != -9999)

    aData_out = aData_clip[dummy_index0]

    print(aData_out)

if __name__ == '__main__':
    sWorkspace_dem_normal = '/compyfs/liao313/00raw/dem/global/gebco/normal'
    sWorkspace_dem_subice = '/compyfs/liao313/00raw/dem/global/gebco/sub_ice'
    aFilename_dem_geotiff = list()
    aFilename_dem_geotiff.append(os.path.join(sWorkspace_dem_normal, 'gebco_2024_n0.0_s-90.0_w0.0_e90.0.tif'))
    aFilename_dem_geotiff.append(os.path.join(sWorkspace_dem_normal, 'gebco_2024_n0.0_s-90.0_w-180.0_e-90.0.tif'))
    aFilename_dem_geotiff.append(os.path.join(sWorkspace_dem_normal, 'gebco_2024_n0.0_s-90.0_w-90.0_e0.0.tif'))
    aFilename_dem_geotiff.append(os.path.join(sWorkspace_dem_normal, 'gebco_2024_n0.0_s-90.0_w90.0_e180.0.tif'))
    aFilename_dem_geotiff.append(os.path.join(sWorkspace_dem_normal, 'gebco_2024_n90.0_s0.0_w0.0_e90.0.tif'))
    aFilename_dem_geotiff.append(os.path.join(sWorkspace_dem_normal, 'gebco_2024_n90.0_s0.0_w-180.0_e-90.0.tif'))
    aFilename_dem_geotiff.append(os.path.join(sWorkspace_dem_normal, 'gebco_2024_n90.0_s0.0_w-90.0_e0.0.tif'))
    aFilename_dem_geotiff.append(os.path.join(sWorkspace_dem_normal, 'gebco_2024_n90.0_s0.0_w90.0_e180.0.tif'))

    sFilename_cut = '/qfs/people/liao313/workspace/python/hexwatershed_utility/figures/mpas/cross_antimeridian.shp'

    gdal_warp_antimeridian_polygon(sFilename_cut, aFilename_dem_geotiff)
