#this function injects the DEM into the mesh file, most of the time, the mesh should be the MPAS mesh
import os, sys
import math
import importlib.util
import contextlib
import numpy as np
import netCDF4 as nc
from osgeo import gdal, osr, ogr
from pyearth.gis.gdal.read.raster.gdal_read_geotiff_file import gdal_read_geotiff_file
from pyearth.toolbox.analysis.extract.clip_raster_by_polygon_file import clip_raster_by_polygon_file
from pyearth.gis.gdal.gdal_validate_polygon_file import gdal_validate_polygon_file
gdal.UseExceptions()
gdal.PushErrorHandler('CPLQuietErrorHandler')
iFlag_cython = importlib.util.find_spec("cython")
if iFlag_cython is not None:
    from tinyr import RTree

from pyearth.gis.geometry.convert_longitude_range import convert_360_to_180_np
from pyearth.toolbox.data.geoparquet.convert_geojson_to_geoparquet import convert_geojson_to_geoparquet
def convert_mpas_mesh_to_geojson( sFilename_mpas_mesh_netcdf_cull,
                                sFilename_geojson_out, sFilename_error_geojson_out = None,
                                sFilename_mpas_mesh_netcdf_base = None ):

    pDriver_geojson = ogr.GetDriverByName('GeoJSON')
    #check file existence
    if not os.path.exists(sFilename_mpas_mesh_netcdf_cull):
        print("The mesh file does not exist, please check the file path")
        return

    #check output file, remove if exists

    if os.path.exists(sFilename_geojson_out):
        os.remove(sFilename_geojson_out)

    if sFilename_error_geojson_out is not None:
        if os.path.exists(sFilename_error_geojson_out):
            os.remove(sFilename_error_geojson_out)

    #first, we need to read the mesh file and copy it to the new file
    pDatasets_mesh = nc.Dataset(sFilename_mpas_mesh_netcdf_cull, 'r')
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
    #aVertexOnEdge0 = verticesOnEdge0[:]
    aIndexToCellID = indexToCellID0[:]
    ncell = len(aIndexToCellID)

    if sFilename_mpas_mesh_netcdf_base is not None:
        pDatasets_mesh_base = nc.Dataset(sFilename_mpas_mesh_netcdf_base, 'r')
        for sKey, aValue in pDatasets_mesh_base.variables.items():
            if sKey == 'lonCell':
                lonCell1 = aValue

            if sKey == 'latCell':
                latCell1 = aValue

            if sKey == 'verticesOnCell':
                verticesOnCell1 = aValue

            if sKey == 'indexToCellID':
                indexToCellID1 = aValue

            if sKey == 'lonVertex':
                lonVertex1 = aValue

            if sKey == 'latVertex':
                latVertex1 = aValue

        aLongitudeCell_base = lonCell1[:] / math.pi * 180
        aLatitudeCell_base = latCell1[:] / math.pi * 180
        aLongitudeCell_base = np.array(aLongitudeCell_base)
        aLatitudeCell_base = np.array(aLatitudeCell_base)
        aLongitudeVertex_base = lonVertex1[:] / math.pi * 180
        aLatitudeVertex_base = latVertex1[:] / math.pi * 180
        aVertexOnCell_base = verticesOnCell1[:]
        aIndexToCellID_base = indexToCellID1[:]
        ncell_base = len(aIndexToCellID_base)
        aIndexToCellID_base = np.array(aIndexToCellID_base)
        aLongitudeCell_base_180 = convert_360_to_180_np(aLongitudeCell_base)
        aLongitudeVertex_base_180 = convert_360_to_180_np(aLongitudeVertex_base)

        index_vertex = RTree(max_cap=5, min_cap=2)
        for i in range(ncell_base):
            lCellID = int(aIndexToCellID_base[i])
            x = float(aLongitudeCell_base_180[i])
            y = float(aLatitudeCell_base[i])
            left =   x - 1E-5
            right =  x + 1E-5
            bottom = y - 1E-5
            top =    y + 1E-5
            pBound= (left, bottom, right, top)
            index_vertex.insert(lCellID, pBound) #
            pass

    pSpatialRef_target = osr.SpatialReference()
    pSpatialRef_target.ImportFromEPSG(4326)

    pDataset = pDriver_geojson.CreateDataSource(sFilename_geojson_out)
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

    if sFilename_error_geojson_out is not None:
        #create a geojson file for the error mesh cell
        pDataset2 = pDriver_geojson.CreateDataSource(sFilename_error_geojson_out)
        pLayerOut2 = pDataset2.CreateLayer('cell', pSpatialRef_target, ogr.wkbPoint)
        pFieldDefn2 = ogr.FieldDefn('id', ogr.OFTInteger)
        pLayerOut2.CreateField(pFieldDefn2)
        pFieldDefn2 = ogr.FieldDefn('error', ogr.OFTInteger)
        pLayerOut2.CreateField(pFieldDefn2)
        pFieldDefn2 = ogr.FieldDefn('lon', ogr.OFTReal)
        pFieldDefn2.SetWidth(25)
        pFieldDefn2.SetPrecision(15)
        pLayerOut2.CreateField(pFieldDefn2)
        pFieldDefn2 = ogr.FieldDefn('lat', ogr.OFTReal)
        pFieldDefn2.SetWidth(25)
        pFieldDefn2.SetPrecision(15)
        pLayerOut2.CreateField(pFieldDefn2)
        pLayerDefn2 = pLayerOut2.GetLayerDefn()

    pPolygon_empty = ogr.Geometry(ogr.wkbPolygon)

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
            #print(aLonVertex, ' and ', aLatVertex)
            if sFilename_error_geojson_out is not None:
                #create a geojson file for the error mesh cell
                pPoint = ogr.Geometry(ogr.wkbPoint)
                pPoint.AddPoint(dLon_center, dLat_center)
                pFeatureOut2 = ogr.Feature(pLayerDefn2)
                pFeatureOut2.SetGeometry(pPoint)
                pFeatureOut2.SetField('id', lCellID)
                pFeatureOut2.SetField('error', 1)
                pFeatureOut2.SetField('lon', dLon_center)
                pFeatureOut2.SetField('lat', dLat_center)
                pLayerOut2.CreateFeature(pFeatureOut2)
            if sFilename_mpas_mesh_netcdf_base is not None:
                aIntersect = list(index_vertex.search_surrounding([dLon_center, dLat_center]))
                #if len(lCell_index_base) == 1 and lCell_index_base[0] != 0:
                if len(aIntersect) == 1:
                    lCell_index_base = aIntersect[0]-1
                    lCellID_base = aIndexToCellID_base[lCell_index_base]
                    print("find in base mesh: ", dLon_center, dLat_center, lCellID, lCellID_base, aIntersect[0])
                    aVertexOnCellIndex = np.array(aVertexOnCell_base[lCell_index_base,:])
                    dummy0 = np.where(aVertexOnCellIndex > 0)
                    aVertexIndex = aVertexOnCellIndex[dummy0] - 1
                    aLonVertex_base = aLongitudeVertex_base_180[aVertexIndex]
                    aLatVertex_base = aLatitudeVertex_base[aVertexIndex]
                    nVertex = len(aLonVertex_base)
                    ring = ogr.Geometry(ogr.wkbLinearRing)
                    for lon, lat in zip(aLonVertex_base, aLatVertex_base):
                        ring.AddPoint(lon, float(lat))

                    ring.CloseRings()
                    pPolygon = ogr.Geometry(ogr.wkbPolygon)
                    pPolygon.AddGeometry(ring)
                    pFeatureOut = ogr.Feature(pLayerDefn)
                    if not pPolygon.IsValid():
                        print("Polygon WKT:", pPolygon.ExportToWkt())
                        #continue
                        pFeatureOut.SetGeometry(pPolygon_empty)  # No geometry
                    else:
                        pFeatureOut.SetGeometry(pPolygon)
                        pass

                    pFeatureOut.SetField('id', lCellID)
                    pFeatureOut.SetField('lon', dLon_center)
                    pFeatureOut.SetField('lat', dLat_center)
                    pLayerOut.CreateFeature(pFeatureOut)
                else:
                    lCell_index_base0 = np.where(aLongitudeCell_base_180 == dLon_center)
                    lCell_index_base1 = np.where(aLatitudeCell_base == dLat_center)
                    print("Cannot find the base cell for the error cell type 1: ", dLon_center, dLat_center, lCellID)
                    print(aLatitudeCell_base[lCell_index_base0])
                    print(aLongitudeCell_base_180[lCell_index_base1])
        else:
            ring = ogr.Geometry(ogr.wkbLinearRing)
            for lon, lat in zip(aLonVertex, aLatVertex):
                ring.AddPoint(float(lon), float(lat))

            ring.CloseRings()
            pPolygon = ogr.Geometry(ogr.wkbPolygon)
            pPolygon.AddGeometry(ring)
            # Validate the geometry
            if not pPolygon.IsValid():
                #print("Polygon is invalid...", i)
                #print("Polygon WKT:", pPolygon.ExportToWkt())
                if sFilename_error_geojson_out is not None:
                    #create a geojson file for the error mesh cell
                    pPoint = ogr.Geometry(ogr.wkbPoint)
                    pPoint.AddPoint(dLon_center, dLat_center)
                    pFeatureOut2 = ogr.Feature(pLayerDefn2)
                    pFeatureOut2.SetGeometry(pPoint)
                    pFeatureOut2.SetField('id', lCellID)
                    pFeatureOut2.SetField('error', 2)
                    pFeatureOut2.SetField('lon', dLon_center)
                    pFeatureOut2.SetField('lat', dLat_center)
                    pLayerOut2.CreateFeature(pFeatureOut2)

                if sFilename_mpas_mesh_netcdf_base is not None:
                    aIntersect = list(index_vertex.search_surrounding([dLon_center, dLat_center]))
                    #find the base mesh cell using the center lon and lat
                    #lCell_index_base = np.where((aLongitudeCell_base_180 == dLon_center) & (aLatitudeCell_base == dLat_center))
                    #check we can find a base cell or not
                    pFeatureOut = ogr.Feature(pLayerDefn)
                    if len(aIntersect) == 1: #and lCell_index_base[0] != 0:
                        lCell_index_base = aIntersect[0]-1
                        lCellID_base = aIndexToCellID_base[lCell_index_base]
                        aVertexOnCellIndex = np.array(aVertexOnCell_base[lCell_index_base,:])
                        dummy0 = np.where(aVertexOnCellIndex > 0)
                        aVertexIndex = aVertexOnCellIndex[dummy0] - 1
                        aLonVertex_base = aLongitudeVertex_base_180[aVertexIndex]
                        aLatVertex_base = aLatitudeVertex_base[aVertexIndex]
                        nVertex = len(aLonVertex_base)
                        ring = ogr.Geometry(ogr.wkbLinearRing)
                        for lon, lat in zip(aLonVertex_base, aLatVertex_base):
                            ring.AddPoint(lon, float(lat))

                        ring.CloseRings()
                        pPolygon = ogr.Geometry(ogr.wkbPolygon)
                        pPolygon.AddGeometry(ring)
                        if not pPolygon.IsValid():
                            print("Polygon WKT:", pPolygon.ExportToWkt())
                            pFeatureOut.SetGeometry(pPolygon_empty)
                            #continue
                        else:
                            pFeatureOut.SetGeometry(pPolygon)
                    else:
                        #continue
                        pFeatureOut.SetGeometry(pPolygon_empty)  # No geometry
                        pass

                    pFeatureOut.SetField('id', lCellID)
                    pFeatureOut.SetField('lon', dLon_center)
                    pFeatureOut.SetField('lat', dLat_center)
                    pLayerOut.CreateFeature(pFeatureOut)

            else:
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
    if sFilename_error_geojson_out is not None:
        pDataset2.FlushCache()
        pFeatureOut2 = None
        pLayerOut2 = None
        pDataset2 = None

    #check whether all the polygons are valid
    iFlag_valid = gdal_validate_polygon_file(sFilename_geojson_out)
    print("The validity of the polygon is: ", iFlag_valid)

    return

if __name__ == '__main__':
    sFilename_mpas_mesh_netcdf_cull = '/compyfs/liao313/04model/pyflowline/global/pyflowline20250101007/jigsaw/out/invert_mesh.nc'
    sFilename_mpas_mesh_netcdf_base = '/compyfs/liao313/04model/pyflowline/global/pyflowline20250101007/jigsaw/out/base_mesh.nc'

    sFilename_geojson_out = '/compyfs/liao313/00raw/mesh/global/hrm/mpas_land_with_empty_geometry.geojson'
    sFilename_error_geojson_out = '/compyfs/liao313/00raw/mesh/global/hrm/mpas_land_error.geojson'

    convert_mpas_mesh_to_geojson(  sFilename_mpas_mesh_netcdf_cull,
                                sFilename_geojson_out,
                                 sFilename_error_geojson_out = sFilename_error_geojson_out,
                                 sFilename_mpas_mesh_netcdf_base = sFilename_mpas_mesh_netcdf_base )

    #convert to geoparquet
    sFilename_geoparquet_out = '/compyfs/liao313/00raw/mesh/global/hrm/mpas_land.parquet'
    convert_geojson_to_geoparquet(sFilename_geojson_out, sFilename_geoparquet_out)