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
if iFlag_cython is not None:
    from pyflowline.algorithms.cython.kernel import convert_360_to_180
else:
    from pyearth.gis.geometry.convert_longitude_range import convert_360_to_180
gdal.UseExceptions()

from hexwatershed_utility.mpas.map_wkt_polygon import map_wkt_polygon

def check_mesh_cell_polygon(sFilename_mpas_mesh_netcdf_culled, aCellID_in,
                             iFlag_plot_in = None ,
                            sFilename_mpas_mesh_netcdf_base = None):
    #sFilename_mpas_mesh_netcdf_culled, the culled mesh cell netcdf
    #aCellID_in, the list of cell ID to be checked
    #iFlag_plot_in, if it is not None, it will plot the cell polygon
    #sFilename_mpas_mesh_netcdf_base, the base mesh cell netcdf, it will check whether the base cell is valid or not, using mesh cell center
    #first, we need to read the mesh file and copy it to the new file
    pDatasets_mesh = nc.Dataset(sFilename_mpas_mesh_netcdf_culled, 'r')
    #get netcdf format
    format = pDatasets_mesh.file_format
    print('Netcdf format:', format)
    #read new netcdf
    for sKey, aValue in pDatasets_mesh.variables.items():
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
    aLongitudeVertex = lonVertex0[:] / math.pi * 180
    aLatitudeVertex = latVertex0[:] / math.pi * 180
    aVertexOnCell = verticesOnCell0[:]
    aIndexToCellID = indexToCellID0[:]
    ncell = len(aIndexToCellID)
    aIndexToCellID = np.array(aIndexToCellID)

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


    for lCellID_in in aCellID_in:
        if lCellID_in > ncell or lCellID_in < 0:
            print('Error: the cell ID is out of range')
            return

        #use the cell ID to obtain its index, all MPAS index starts from 1!

        lCellIndex = np.where(aIndexToCellID == lCellID_in)
        dLon_center_360 = float(aLongitudeCell[lCellIndex])
        dLon_center = convert_360_to_180 (dLon_center_360)
        dLat_center =  float(aLatitudeCell[lCellIndex])
        aVertexOnCellIndex = np.array(aVertexOnCell[lCellIndex,:])
        dummy0 = np.where(aVertexOnCellIndex > 0)
        aVertexIndex = aVertexOnCellIndex[dummy0] - 1
        aLonVertex = aLongitudeVertex[aVertexIndex]
        aLatVertex = aLatitudeVertex[aVertexIndex]
        nVertex = len(aLonVertex)
        if nVertex < 3:
            print("Vertex number is: ", nVertex)
            return
        else:
            ring = ogr.Geometry(ogr.wkbLinearRing)
            aCoords_gcs = np.full((nVertex,2), -9999.0, dtype=float)
            for j in range(nVertex):
                x1 = convert_360_to_180(aLonVertex[j])
                y1 = float(aLatVertex[j])
                ring.AddPoint(x1, y1)
                aCoords_gcs[j,0] = x1
                aCoords_gcs[j,1] = y1
                pass

            #add the closing point, which is the first vertex
            x1 = convert_360_to_180(aLonVertex[0])
            y1 = aLatVertex[0]
            ring.AddPoint(x1, y1) #double check
            pPolygon = ogr.Geometry(ogr.wkbPolygon)
            pPolygon.AddGeometry(ring)
            # Validate the geometry
            pPolygon_wkt = pPolygon.ExportToWkt()
            if not pPolygon.IsValid():
                print("Polygon is invalid...")
                #pPolygon = pPolygon.Buffer(0)  # Buffering by 0 can fix some geometry issues
                #pPolygon = pPolygon.Simplify(0.0001)  # Simplify with a tolerance
            else:
                print("Polygon is valid...")

            print("Polygon WKT:", pPolygon_wkt)
            if iFlag_plot_in is not None:
                #use a simple method to plot the polygon
                #get the path where this script is located
                sPath_script = os.path.dirname(os.path.abspath(__file__))
                sFilename_png = 'mesh_cell_' + str(lCellID_in) + '_culled.png'
                sFilename_png_out = os.path.join(sPath_script, sFilename_png)
                map_wkt_polygon(pPolygon_wkt, sFilename_png_out)
                if sFilename_mpas_mesh_netcdf_base is not None:
                    #find the base mesh cell using the center lon and lat
                    lCell_index_base = np.where((aLongitudeCell_base == dLon_center_360) & (aLatitudeCell_base == dLat_center))
                    #check we can find a base cell or not
                    if len(lCell_index_base) == 1 and lCell_index_base[0] != 0:
                        lCellID_base = aIndexToCellID_base[lCell_index_base]
                        aVertexOnCellIndex = np.array(aVertexOnCell_base[lCell_index_base,:])
                        dummy0 = np.where(aVertexOnCellIndex > 0)
                        aVertexIndex = aVertexOnCellIndex[dummy0] - 1
                        aLonVertex_base = aLongitudeVertex_base[aVertexIndex]
                        aLatVertex_base = aLatitudeVertex_base[aVertexIndex]
                        nVertex = len(aLonVertex)
                        if nVertex < 3:
                            print("Base vertex number is: ", nVertex)
                            return
                        else:
                            ring = ogr.Geometry(ogr.wkbLinearRing)
                            aCoords_gcs = np.full((nVertex,2), -9999.0, dtype=float)
                            for j in range(nVertex):
                                x1 = convert_360_to_180(aLonVertex_base[j])
                                y1 = float(aLatVertex_base[j])
                                ring.AddPoint(x1, y1)
                                aCoords_gcs[j,0] = x1
                                aCoords_gcs[j,1] = y1
                                pass

                            #add the closing point, which is the first vertex
                            x1 = convert_360_to_180(aLonVertex_base[0])
                            y1 = aLatVertex_base[0]
                            ring.AddPoint(x1, y1) #double check
                            pPolygon = ogr.Geometry(ogr.wkbPolygon)
                            pPolygon.AddGeometry(ring)
                            # Validate the geometry
                            pPolygon_wkt_base = pPolygon.ExportToWkt()
                            if not pPolygon.IsValid():
                                print("Base polygon is invalid...")
                                #pPolygon = pPolygon.Buffer(0)  # Buffering by 0 can fix some geometry issues
                                #pPolygon = pPolygon.Simplify(0.0001)  # Simplify with a tolerance
                            else:
                                print("Base polygon is valid...")

                            sFilename_png = 'mesh_cell_' + str(lCellID_in) + '_base.png'
                            sFilename_png_out = os.path.join(sPath_script, sFilename_png)
                            map_wkt_polygon(pPolygon_wkt_base, sFilename_png_out)
                    else:
                        print("The base cell is not found...")
                        continue
                pass



    return

#create a main call
if __name__ == '__main__':
    sFilename_mpas_mesh_netcdf_culled = '/compyfs/liao313/04model/pyflowline/conus/pyflowline20241201019/jigsaw/out/invert_mesh.nc'
    sFilename_mpas_mesh_netcdf_base = '/compyfs/liao313/04model/pyflowline/conus/pyflowline20241201019/jigsaw/out/base_mesh.nc'
   


    check_mesh_cell_polygon(sFilename_mpas_mesh_netcdf_culled, aCellID_in,  iFlag_plot_in = 1, sFilename_mpas_mesh_netcdf_base = sFilename_mpas_mesh_netcdf_base )