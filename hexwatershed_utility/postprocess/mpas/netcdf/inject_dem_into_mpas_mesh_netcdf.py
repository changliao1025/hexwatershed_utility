#this function injects the DEM into the mesh file, most of the time, the mesh should be the MPAS mesh
import os, sys
import math
import importlib.util
import numpy as np
import netCDF4 as nc
from osgeo import gdal, osr, ogr
from pyearth.gis.gdal.read.raster.gdal_read_geotiff_file import gdal_read_geotiff_file
from pyearth.toolbox.analysis.extract.clip_raster_by_polygon_file import clip_raster_by_polygon_file
from pyearth.gis.geometry.convert_longitude_range import convert_360_to_180_np
from pyearth.gis.geometry.convert_idl_polygon_to_valid_polygon import convert_idl_polygon_to_valid_polygon
from pyearth.gis.geometry.split_polygon_cross_idl import split_polygon_cross_idl

gdal.UseExceptions()

def inject_dem_into_mpas_mesh_netcdf(  sFilename_mesh_netcdf,
                                sFilename_mesh_with_dem_netcdf,
                                aFilename_dem_geotiff,
                                aFilename_dem_with_ice ):

    pDriver_shapefile = ogr.GetDriverByName('ESRI Shapefile')
    #check file existence
    if not os.path.exists(sFilename_mesh_netcdf):
        print("The mesh file does not exist, please check the file path")
        return

    for sFilename_dem_geotiff in aFilename_dem_geotiff:
        if not os.path.exists(sFilename_dem_geotiff):
            print("The DEM file does not exist, please check the file path")
            return


    for sFilename_dem_with_ice in aFilename_dem_with_ice:
        if not os.path.exists(sFilename_dem_with_ice):
            print("The DEM with ice file does not exist, please check the file path")

    #check output file, remove if exists
    if os.path.exists(sFilename_mesh_with_dem_netcdf):
        os.remove(sFilename_mesh_with_dem_netcdf)

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

    aFilename_dem_with_ice_masked = list()
    for sFilename_dem_with_ice in aFilename_dem_with_ice:
        sFilename_dem_with_ice_masked = sFilename_dem_with_ice.replace('.tif', '_masked.tif')
        #if os.path.exists(sFilename_dem_with_ice_masked):
        #    os.remove(sFilename_dem_with_ice_masked)
        #clip_raster_by_polygon_file(sFilename_dem_with_ice, sFilename_land_ocean_mask, sFilename_dem_with_ice_masked,
        #                            iFlag_use_raster_extent = 1 )
        aFilename_dem_with_ice_masked.append(sFilename_dem_with_ice_masked)


    #first, we need to read the mesh file and copy it to the new file
    pDatasets_mesh = nc.Dataset(sFilename_mesh_netcdf, 'r')
    #get netcdf format
    format = pDatasets_mesh.file_format
    #create a new file
    pDatasets_out = nc.Dataset(sFilename_mesh_with_dem_netcdf, 'w',format=format)
    #use a simliar approach in pyearth to copy the mesh file
    aDimension_key=list()
    aDimension_value=list()
    for sKey, iValue in pDatasets_mesh.dimensions.items():
        dummy = len(iValue)
        if not iValue.isunlimited():
            aDimension_key.append(sKey)
            aDimension_value.append(sKey)
            pDatasets_out.createDimension(sKey, dummy)
        else:
            pDatasets_out.createDimension(sKey, dummy )

    #read new netcdf
    for sKey, aValue in pDatasets_mesh.variables.items():
        #we need to filter out unused grids based on mpas specs
        if sKey == 'verticesOnCell':
            verticesOnCell0 = aValue

        if sKey == 'indexToCellID':
            indexToCellID0 = aValue

        if sKey == 'lonVertex':
            lonVertex0 = aValue

        if sKey == 'latVertex':
            latVertex0 = aValue

        if sKey == 'lonCell':
            lonCell0 = aValue

        if sKey == 'latCell':
            latCell0 = aValue

    aLatitudeVertex = latVertex0[:] / math.pi * 180
    aLongitudeVertex = lonVertex0[:] / math.pi * 180
    #convert unit
    aLatitudeCell = latCell0[:] / math.pi * 180
    aLongitudeCell = lonCell0[:] / math.pi * 180
    #aCellsOnCell = cellsOnCell0[:]
    #aCellOnEdge = cellsOnEdge0[:]
    #aEdgesOnCell = edgesOnCell0[:]
    aVertexOnCell = verticesOnCell0[:]
    #aVertexOnEdge0 = verticesOnEdge0[:]
    aIndexToCellID = indexToCellID0[:]
    ncell = len(aIndexToCellID)
    aLongitudeCell_base_180 = convert_360_to_180_np(aLongitudeCell)
    aLongitudeVertex_base_180 = convert_360_to_180_np(aLongitudeVertex)
    #add the new variable dimension?
    #find out which dimension equal to nCells
    for sKey, iValue in pDatasets_mesh.dimensions.items():
        dummy = len(iValue)
        if dummy == ncell:
            sDimension_to_be_add = sKey
            break

    for sKey, aValue in pDatasets_mesh.variables.items():
        # we need to take care of rec dimension
        dummy = aValue.dimensions
        #check is the fill value exist or not
        if '_FillValue' in aValue.ncattrs():
            outVar = pDatasets_out.createVariable(sKey, aValue.datatype, dummy, fill_value=aValue._FillValue)
        else:
            outVar = pDatasets_out.createVariable(sKey, aValue.datatype, dummy)

        for sAttribute in aValue.ncattrs():
            if sAttribute != '_FillValue':
                outVar.setncatts( { sAttribute: aValue.getncattr(sAttribute) } )

        outVar[:] = aValue[:]

    #now we need to read the DEM file, using the subgrid method

    dummy = gdal_read_geotiff_file(sFilename_dem_masked)
    dPixelWidth = dummy['pixelWidth']
    pPixelHeight = dummy['pixelHeight']
    dMissing_value= dummy['missingValue']
    pProjection = dummy['projection']
    pSpatialRef_target= osr.SpatialReference()
    pSpatialRef_target.ImportFromWkt(pProjection)
    dummy = gdal_read_geotiff_file(sFilename_dem_with_ice_masked)
    dMissing_value1= dummy['missingValue']
    #get the DEM information
    sFilename_shapefile_cut = "/vsimem/tmp_polygon.shp"
    sResampleAlg = 'near'
    aData_bed_elevation = np.full((ncell), -9999.0, dtype=float)
    aData_ice_thickness =  np.full((ncell), -9999.0, dtype=float)
    #debug
    iFlag_debug = 0
    if iFlag_debug == 1:
        i_start = 303141
    else:
        i_start = 0
    for i in range(i_start, ncell, 1):
        #center
        dLon = aLongitudeCell[i]
        dLat = aLatitudeCell[i]
        aVertexOnCellIndex = np.array(aVertexOnCell[i,:])
        dummy0 = np.where(aVertexOnCellIndex > 0)
        aVertexIndex = aVertexOnCellIndex[dummy0]
        aLonVertex = aLongitudeVertex_base_180[aVertexIndex-1]
        aLatVertex = aLatitudeVertex[aVertexIndex-1]
        nVertex = len(aLonVertex)
        if nVertex < 3:
            print("Vertex number is: ", nVertex, i)
            continue
        else:
            #ignore antarctica for now
            if dLat < -60:
                continue
            #also ignore the north pole
            if dLat > 85:
                continue

            #first check if it is within the boundary
            ring = ogr.Geometry(ogr.wkbLinearRing)
            aCoords_gcs = np.full((nVertex,2), -9999.0, dtype=float)
            for j in range(nVertex):
                x1 = aLonVertex[j]
                y1 = aLatVertex[j]
                ring.AddPoint(x1, y1)
                aCoords_gcs[j,0] = x1
                aCoords_gcs[j,1] = y1
                pass

            x1 = aLonVertex[0]
            y1 = aLatVertex[0]
            ring.AddPoint(x1, y1) #double check
            pPolygon = ogr.Geometry(ogr.wkbPolygon)
            pPolygon.AddGeometry(ring)
            dLon_min = np.min(aCoords_gcs[:,0])
            dLon_max = np.max(aCoords_gcs[:,0])
            if np.abs(dLon_min-dLon_max) > 100:
                #cross the internation date line
                pPolygon_new = convert_idl_polygon_to_valid_polygon(pPolygon)
                if pPolygon_new is not None:
                    if pPolygon_new.IsValid() == False:
                        print('Warning: invalid polygon')
                        continue
                    else:
                        pPolygon_new.FlattenTo2D()
                        print("Polygon WKT:", pPolygon_new.ExportToWkt())
                        print('Splitting polygon...', i)
                        #can be split and then used for clipping
                        aCoord_gcs_split = split_polygon_cross_idl(aCoords_gcs)
                        #if iFlag_debug == 1:
                        #    sFilename_shapefile_cut = '/qfs/people/liao313/workspace/python/hexwatershed_utility/hexwatershed_utility/mpas/netcdf/debug/tmp_polygon.shp'
                        #    if os.path.exists(sFilename_shapefile_cut):
                        #        pDriver_shapefile.DeleteDataSource(sFilename_shapefile_cut)
                        pDataset3 = pDriver_shapefile.CreateDataSource(sFilename_shapefile_cut)
                        pLayerOut3 = pDataset3.CreateLayer('cell', pSpatialRef_target, ogr.wkbPolygon)
                        pLayerDefn3 = pLayerOut3.GetLayerDefn()
                        for aCoord_gcs in aCoord_gcs_split:
                            ring = ogr.Geometry(ogr.wkbLinearRing)
                            for aCoord in aCoord_gcs:
                                x1 = aCoord[0]
                                y1 = aCoord[1]
                                ring.AddPoint(x1, y1)
                                pass
                            x1 = aCoord_gcs[0][0]
                            y1 = aCoord_gcs[0][1]
                            ring.AddPoint(x1, y1)
                            pPolygon_new = ogr.Geometry(ogr.wkbPolygon)
                            pPolygon_new.AddGeometry(ring)
                            pFeatureOut3 = ogr.Feature(pLayerDefn3)
                            pFeatureOut3.SetGeometry(pPolygon_new)
                            pLayerOut3.CreateFeature(pFeatureOut3)
                        pDataset3.FlushCache()
                        pass
                else:
                    continue
                pass
            else:
                # Validate the geometry
                if not pPolygon.IsValid():
                    print("Polygon is invalid...", i)
                    pPolygon.FlattenTo2D()
                    print("Polygon WKT:", pPolygon.ExportToWkt())
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
                dummy = aData_clip[dummy_index0]
                if dummy.size > 0:
                    aData_bed_elevation[i] = float(np.mean(dummy))
                else:
                    print('bed elevation issue', i)

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
                dummy = aIce_thickness[dummy_index1]
                if dummy.size > 0:
                    aData_ice_thickness[i] = float(np.mean(dummy))
                else:
                    print('ice thickness issue', i)

            #print('clip success', i)


    pVar = pDatasets_out.createVariable('bed_elevation', 'f8', (sDimension_to_be_add,))
    pVar[:] = aData_bed_elevation
    pVar.setncatts( { 'units': 'm' } )
    pVar.setncatts( { 'long_name': 'bed elevation' } )


    aData_ice_thickness = np.array(aData_ice_thickness)
    pVar = pDatasets_out.createVariable('ice_thickness', 'f8', (sDimension_to_be_add,))
    pVar[:] = aData_ice_thickness
    pVar.setncatts( { 'units': 'm' } )
    pVar.setncatts( { 'long_name': 'ice thickness' } )

    #close the dataset
    pDatasets_out.close()
    pDatasets_mesh.close()

    return

if __name__ == '__main__':
    sFilename_mesh_netcdf = '/compyfs/liao313/04model/pyflowline/global/pyflowline20250101010/jigsaw/out/invert_mesh.nc'
    sFilename_mesh_with_dem_netcdf='/compyfs/liao313/00raw/mesh/global/hrm/mpas_with_land_dem3.nc'
    sWorkspace_dem_normal = '/compyfs/liao313/00raw/dem/global/gebco/normal'
    sWorkspace_dem_subice = '/compyfs/liao313/00raw/dem/global/gebco/sub_ice'
    #gebco_2024_n0.0_s-90.0_w0.0_e90.0.tif*
    #gebco_2024_n0.0_s-90.0_w-180.0_e-90.0.tif*
    #gebco_2024_n0.0_s-90.0_w-90.0_e0.0.tif*
    #gebco_2024_n0.0_s-90.0_w90.0_e180.0.tif*
    #gebco_2024_n90.0_s0.0_w0.0_e90.0.tif*
    #gebco_2024_n90.0_s0.0_w-180.0_e-90.0.tif*
    #gebco_2024_n90.0_s0.0_w-90.0_e0.0.tif*
    #gebco_2024_n90.0_s0.0_w90.0_e180.0.tif*

    aFilename_dem_geotiff = list()
    aFilename_dem_geotiff.append(os.path.join(sWorkspace_dem_normal, 'gebco_2024_n0.0_s-90.0_w0.0_e90.0.tif'))
    aFilename_dem_geotiff.append(os.path.join(sWorkspace_dem_normal, 'gebco_2024_n0.0_s-90.0_w-180.0_e-90.0.tif'))
    aFilename_dem_geotiff.append(os.path.join(sWorkspace_dem_normal, 'gebco_2024_n0.0_s-90.0_w-90.0_e0.0.tif'))
    aFilename_dem_geotiff.append(os.path.join(sWorkspace_dem_normal, 'gebco_2024_n0.0_s-90.0_w90.0_e180.0.tif'))
    aFilename_dem_geotiff.append(os.path.join(sWorkspace_dem_normal, 'gebco_2024_n90.0_s0.0_w0.0_e90.0.tif'))
    aFilename_dem_geotiff.append(os.path.join(sWorkspace_dem_normal, 'gebco_2024_n90.0_s0.0_w-180.0_e-90.0.tif'))
    aFilename_dem_geotiff.append(os.path.join(sWorkspace_dem_normal, 'gebco_2024_n90.0_s0.0_w-90.0_e0.0.tif'))
    aFilename_dem_geotiff.append(os.path.join(sWorkspace_dem_normal, 'gebco_2024_n90.0_s0.0_w90.0_e180.0.tif'))

    #gebco_2024_sub_ice_n0.0_s-90.0_w0.0_e90.0.tif*
    #gebco_2024_sub_ice_n0.0_s-90.0_w-180.0_e-90.0.tif*
    #gebco_2024_sub_ice_n0.0_s-90.0_w-90.0_e0.0.tif*
    #gebco_2024_sub_ice_n0.0_s-90.0_w90.0_e180.0.tif*
    #gebco_2024_sub_ice_n90.0_s0.0_w0.0_e90.0.tif*
    #gebco_2024_sub_ice_n90.0_s0.0_w-180.0_e-90.0.tif*
    #gebco_2024_sub_ice_n90.0_s0.0_w-90.0_e0.0.tif*
    #gebco_2024_sub_ice_n90.0_s0.0_w90.0_e180.0.tif*
    aFilename_dem_with_ice = list()
    aFilename_dem_with_ice.append(os.path.join(sWorkspace_dem_subice, 'gebco_2024_sub_ice_n0.0_s-90.0_w0.0_e90.0.tif'))
    aFilename_dem_with_ice.append(os.path.join(sWorkspace_dem_subice, 'gebco_2024_sub_ice_n0.0_s-90.0_w-180.0_e-90.0.tif'))
    aFilename_dem_with_ice.append(os.path.join(sWorkspace_dem_subice, 'gebco_2024_sub_ice_n0.0_s-90.0_w-90.0_e0.0.tif'))
    aFilename_dem_with_ice.append(os.path.join(sWorkspace_dem_subice, 'gebco_2024_sub_ice_n0.0_s-90.0_w90.0_e180.0.tif'))
    aFilename_dem_with_ice.append(os.path.join(sWorkspace_dem_subice, 'gebco_2024_sub_ice_n90.0_s0.0_w0.0_e90.0.tif'))
    aFilename_dem_with_ice.append(os.path.join(sWorkspace_dem_subice, 'gebco_2024_sub_ice_n90.0_s0.0_w-180.0_e-90.0.tif'))
    aFilename_dem_with_ice.append(os.path.join(sWorkspace_dem_subice, 'gebco_2024_sub_ice_n90.0_s0.0_w-90.0_e0.0.tif'))
    aFilename_dem_with_ice.append(os.path.join(sWorkspace_dem_subice, 'gebco_2024_sub_ice_n90.0_s0.0_w90.0_e180.0.tif'))

    inject_dem_into_mpas_mesh_netcdf(  sFilename_mesh_netcdf,
                                sFilename_mesh_with_dem_netcdf,
                                aFilename_dem_geotiff,
                                aFilename_dem_with_ice )