#this function injects the DEM into the mesh file, most of the time, the mesh should be the MPAS mesh
import os, sys
import numpy as np
from osgeo import gdal, osr, ogr
from pyearth.gis.gdal.read.raster.gdal_read_geotiff_file import gdal_read_geotiff_file
gdal.UseExceptions()

def inject_dem_into_mpas_mesh_geojson( sFilename_geojson_in,
                                 sFilename_geojson_out,
                                aFilename_dem_geotiff,
                                aFilename_dem_with_ice = None ):

    pDriver_geojson = ogr.GetDriverByName('GeoJSON')
    pDriver_shapefile = ogr.GetDriverByName('ESRI Shapefile')
    #check file existence
    if not os.path.exists(sFilename_geojson_in):
        print("The mesh file does not exist, please check the file path")
        return

    if os.path.exists(sFilename_geojson_out):
        os.remove(sFilename_geojson_out)

    for sFilename_dem_geotiff in aFilename_dem_geotiff:
        if not os.path.exists(sFilename_dem_geotiff):
            print("The DEM file does not exist, please check the file path")
            return

    if aFilename_dem_with_ice is not None:
        for sFilename_dem_with_ice in aFilename_dem_with_ice:
            if not os.path.exists(sFilename_dem_with_ice):
                print("The DEM with ice file does not exist, please check the file path")

    #first we need to extract the DEM using land ocean mask
    sFilename_land_ocean_mask = '/qfs/people/liao313/data/hexwatershed/global/vector/land_ocean_mask.geojson'

    aFilename_dem_masked = list()
    for sFilename_dem_geotiff in aFilename_dem_geotiff:
        sFilename_dem_masked = sFilename_dem_geotiff.replace('.tif', '_masked.tif')
        #if os.path.exists(sFilename_dem_masked):
        #    os.remove(sFilename_dem_masked)
        #clip_raster_by_polygon_file(sFilename_dem_geotiff, sFilename_land_ocean_mask, sFilename_dem_masked,
        #                            iFlag_use_raster_extent =1)
        aFilename_dem_masked.append(sFilename_dem_masked)
    if aFilename_dem_with_ice is not None:
        aFilename_dem_with_ice_masked = list()
        for sFilename_dem_with_ice in aFilename_dem_with_ice:
            sFilename_dem_with_ice_masked = sFilename_dem_with_ice.replace('.tif', '_masked.tif')
            #if os.path.exists(sFilename_dem_with_ice_masked):
            #    os.remove(sFilename_dem_with_ice_masked)
            #clip_raster_by_polygon_file(sFilename_dem_with_ice, sFilename_land_ocean_mask, sFilename_dem_with_ice_masked,
            #                            iFlag_use_raster_extent = 1 )
            aFilename_dem_with_ice_masked.append(sFilename_dem_with_ice_masked)

    dummy = gdal_read_geotiff_file(sFilename_dem_masked)
    dPixelWidth = dummy['pixelWidth']
    pPixelHeight = dummy['pixelHeight']
    dMissing_value= dummy['missingValue']
    pProjection = dummy['projection']
    pSpatialRef_target= osr.SpatialReference()
    pSpatialRef_target.ImportFromWkt(pProjection)
    dummy = gdal_read_geotiff_file(sFilename_dem_with_ice_masked)
    dMissing_value1= dummy['missingValue']
    sFilename_shapefile_cut = "/vsimem/tmp_polygon.shp"
    sResampleAlg = 'near'

    #create the output geojson file
    pDataset2 = pDriver_geojson.CreateDataSource(sFilename_geojson_out)
    pLayerOut2 = pDataset2.CreateLayer('cell', pSpatialRef_target, ogr.wkbPolygon)

    #add id, area and mean, min, max, std of the raster
    pLayerOut2.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
    pLayerOut2.CreateField(ogr.FieldDefn('bed_elevation', ogr.OFTReal))
    pLayerOut2.CreateField(ogr.FieldDefn('ice_thickness', ogr.OFTReal))
    pLayerDefn2 = pLayerOut2.GetLayerDefn()

    pDataset_subset = ogr.Open(sFilename_geojson_in)
    pLayer_subset = pDataset_subset.GetLayer(0)
    pLayer_subset.ResetReading()
    i = 0
    pPolygon_empty = ogr.Geometry(ogr.wkbPolygon)
    pFeature_subset = pLayer_subset.GetNextFeature()
    while pFeature_subset is not None:
        #use cell id as the name of the output file
        #get cellid from the feature
        lCellId = pFeature_subset.GetField('id')
        sClip = "{:06d}".format(lCellId)
        pFeatureOut2 = ogr.Feature(pLayerDefn2)
        pPolygon = pFeature_subset.GetGeometryRef()
        if pPolygon is None or pPolygon.IsEmpty() or not pPolygon.IsValid():
            print("The polygon is empty")
            #flush stdout
            sys.stdout.flush()
            pFeatureOut2.SetGeometry(pPolygon_empty)
            pFeatureOut2.SetField('id', lCellId)
            pFeatureOut2.SetField('bed_elevation', -9999)
            if aFilename_dem_with_ice is not None:
                pFeatureOut2.SetField('ice_thickness', -9999)
            pLayerOut2.CreateFeature(pFeatureOut2)
            pFeature_subset = pLayer_subset.GetNextFeature()
            i = i + 1
            pDataset2.FlushCache()
            continue
        else:
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
            pDataset_clip_warped = gdal.Warp('', aFilename_dem_masked, options=pWrapOption)
            aData_clip = pDataset_clip_warped.ReadAsArray()
            aData_clip[aData_clip == dMissing_value] = -9999
            dummy_index0 = np.where(aData_clip != -9999)
            if len(dummy_index0) > 0:
                dummy1 = aData_clip[dummy_index0]
                if len(dummy1) > 0:
                    dBed_elevation = float(np.mean(dummy1))
                else:
                    dBed_elevation = -9999
            else:
                dBed_elevation = -9999
            pFeatureOut2.SetGeometry(pPolygon)
            pFeatureOut2.SetField('id', lCellId)
            pFeatureOut2.SetField('bed_elevation', dBed_elevation)
            if aFilename_dem_with_ice is not None:
                pDataset_clip_subice_warped = gdal.Warp('', aFilename_dem_with_ice_masked, options=pWrapOption)
                aData_clip_ice = pDataset_clip_subice_warped.ReadAsArray()
                nan_index0 = np.where(aData_clip == dMissing_value)
                nan_index = np.where(aData_clip_ice == dMissing_value1)
                aIce_thickness =  aData_clip - aData_clip_ice
                aIce_thickness[nan_index0] = -9999
                aIce_thickness[nan_index] = -9999
                dummy_index = np.where(aIce_thickness < 0)
                if len(dummy_index) > 0:
                    aIce_thickness[dummy_index] = 0
                aIce_thickness[nan_index0] = -9999
                aIce_thickness[nan_index] = -9999
                dummy_index1 = np.where(aIce_thickness != -9999)
                if len(dummy_index1) > 0:
                    dummy2 = aIce_thickness[dummy_index1]
                    if len(dummy2) > 0:
                        dIce_thickness = float(np.mean(dummy2))
                    else:
                        dIce_thickness = -9999
                else:
                    dIce_thickness = -9999

                pFeatureOut2.SetField('ice_thickness', dIce_thickness)

            pLayerOut2.CreateFeature(pFeatureOut2)
            pFeature_subset = pLayer_subset.GetNextFeature()
            i = i + 1
            pDataset2.FlushCache()

    #close the dataset
    pDataset2.Destroy()
    pDataset_subset.Destroy()
    pDataset3.Destroy()
    pDataset_clip_warped = None
    print("done")

    return

if __name__ == '__main__':
    sFilename_mesh_geojson_in = '/compyfs/liao313/00raw/mesh/global/hrm/mpas_land_with_empty_geometry.geojson'
    sFilename_mesh_geojson_out='/compyfs/liao313/00raw/mesh/global/hrm/mpas_land_with_elevation.geojson'

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
    aFilename_dem_with_ice = list()
    aFilename_dem_with_ice.append(os.path.join(sWorkspace_dem_subice, 'gebco_2024_sub_ice_n0.0_s-90.0_w0.0_e90.0.tif'))
    aFilename_dem_with_ice.append(os.path.join(sWorkspace_dem_subice, 'gebco_2024_sub_ice_n0.0_s-90.0_w-180.0_e-90.0.tif'))
    aFilename_dem_with_ice.append(os.path.join(sWorkspace_dem_subice, 'gebco_2024_sub_ice_n0.0_s-90.0_w-90.0_e0.0.tif'))
    aFilename_dem_with_ice.append(os.path.join(sWorkspace_dem_subice, 'gebco_2024_sub_ice_n0.0_s-90.0_w90.0_e180.0.tif'))
    aFilename_dem_with_ice.append(os.path.join(sWorkspace_dem_subice, 'gebco_2024_sub_ice_n90.0_s0.0_w0.0_e90.0.tif'))
    aFilename_dem_with_ice.append(os.path.join(sWorkspace_dem_subice, 'gebco_2024_sub_ice_n90.0_s0.0_w-180.0_e-90.0.tif'))
    aFilename_dem_with_ice.append(os.path.join(sWorkspace_dem_subice, 'gebco_2024_sub_ice_n90.0_s0.0_w-90.0_e0.0.tif'))
    aFilename_dem_with_ice.append(os.path.join(sWorkspace_dem_subice, 'gebco_2024_sub_ice_n90.0_s0.0_w90.0_e180.0.tif'))

    inject_dem_into_mpas_mesh_geojson(sFilename_mesh_geojson_in,
                                sFilename_mesh_geojson_out,
                                aFilename_dem_geotiff,
                                aFilename_dem_with_ice = aFilename_dem_with_ice )