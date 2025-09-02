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
from pyearth.gis.geometry.convert_longitude_range import convert_360_to_180_np
from hexwatershed_utility.mpas.map_wkt_polygon import map_wkt_polygon
gdal.UseExceptions()
def exam_mpas_base_mesh(sFilename_mpas_mesh_netcdf_base ,
                             iFlag_plot_in = None , sFolder_out = None):

    if iFlag_plot_in is not None:
        if not os.path.exists(sFolder_out):
            os.makedirs(sFolder_out)

    pDatasets_mesh = nc.Dataset(sFilename_mpas_mesh_netcdf_base, 'r')
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

        if sKey == 'verticesOnEdge':
            verticesOnEdge0 = aValue

        if sKey == 'indexToCellID':
            indexToCellID0 = aValue

        if sKey == 'indexToEdgeID':
            indexToEdgeID0 = aValue

        if sKey == 'lonVertex':
            lonVertex0 = aValue

        if sKey == 'latVertex':
            latVertex0 = aValue

        if sKey == 'edgesOnCell':
            edgesOnCell0 = aValue

    aLongitudeCell = lonCell0[:] / math.pi * 180
    aLatitudeCell = latCell0[:] / math.pi * 180
    aLongitudeVertex = lonVertex0[:] / math.pi * 180
    aLatitudeVertex = latVertex0[:] / math.pi * 180
    aVertexOnCell = verticesOnCell0[:]
    aIndexToCellID = indexToCellID0[:]
    aEdgeOnCell = edgesOnCell0[:]
    aIndexToEdgeID = indexToEdgeID0[:]
    aVerticesOnEdge = verticesOnEdge0[:]

    nCell = len(aIndexToCellID)
    aIndexToCellID = np.array(aIndexToCellID)
    aLongitudeCell_180 = convert_360_to_180_np(aLongitudeCell)
    aLongitudeVertex_180 = convert_360_to_180_np(aLongitudeVertex)

    lCellIndex_debug = 10784
    for i in range(nCell): #lCellIndex_debug, lCellIndex_debug+1, 1):
        #use the cell ID to obtain its index, all MPAS index starts from 1!
        #lCellID_in = i + 1
        lCellIndex = i #np.where(aIndexToCellID == lCellID_in)
        lCellID_in = aIndexToCellID[lCellIndex]
        dLongitude_center =  float(aLongitudeCell_180[i])
        dLatitude_center =  float(aLatitudeCell[i])
        dLat_center =  float(aLatitudeCell[lCellIndex])
        aVertexOnCellIndex = np.array(aVertexOnCell[lCellIndex,:])
        dummy0 = np.where(aVertexOnCellIndex > 0)
        aVertexIndex = aVertexOnCellIndex[dummy0]
        aLonVertex = aLongitudeVertex_180[aVertexIndex-1]
        aLatVertex = aLatitudeVertex[aVertexIndex-1]
        nVertex = len(aLonVertex)
        if nVertex < 3:
            print("Vertex number is: ", nVertex)
            return
        else:
            ring = ogr.Geometry(ogr.wkbLinearRing)
            aCoords_gcs = np.full((nVertex,2), -9999.0, dtype=float)
            for j in range(nVertex):
                x1 = aLonVertex[j]
                y1 = aLatVertex[j]
                ring.AddPoint(x1, y1)
                aCoords_gcs[j,0] = x1
                aCoords_gcs[j,1] = y1
                pass

            #add the closing point, which is the first vertex
            x1 = aLonVertex[0]
            y1 = aLatVertex[0]
            ring.AddPoint(x1, y1) #double check
            pPolygon = ogr.Geometry(ogr.wkbPolygon)
            pPolygon.AddGeometry(ring)
            # Validate the geometry
            #check if the polygon cross the -180 and 180 degree line
            if np.min(aLonVertex) < -150 and np.max(aLonVertex) > 150:
                #pPolygon_wkt = pPolygon.ExportToWkt()
                #print("Polygon crosses the -180 and 180 degree line...", lCellID_in)
                #print("Polygon WKT:", pPolygon_wkt)
                continue

            if not pPolygon.IsValid():
                #157.03394,-7.29670
                #157.03398,-7.29676
                pPolygon.FlattenTo2D()
                pPolygon_wkt = pPolygon.ExportToWkt()
                print("Polygon is invalid...", lCellID_in)
                print("Polygon WKT:", pPolygon_wkt)
                if iFlag_plot_in is not None:
                    #use a simple method to plot the polygon
                    #get the path where this script is located
                    #sPath_script = os.path.dirname(os.path.abspath(__file__))
                    sFilename_png = 'base_mesh_cell_' + str(lCellID_in) + '.png'
                    sFilename_png_out = os.path.join(sFolder_out, sFilename_png)
                    map_wkt_polygon(pPolygon_wkt, sFilename_png_out)

                #check other variables
                #edge on cell
                aEdgeIDOnCell = aEdgeOnCell[lCellIndex,:]
                dummy0 = np.where(aEdgeIDOnCell > 0)
                aEdgeID = aEdgeIDOnCell[dummy0]

                #vertices on edge
                iFlag_print_edge = 0
                if iFlag_print_edge == 1:
                    nEdge = len(aEdgeID)
                    for iEdge in range(nEdge):
                        lEdgeindex = aEdgeID[iEdge] - 1
                        aVertexOnEdgeIndex = aVerticesOnEdge[lEdgeindex,:]
                        print(aVertexOnEdgeIndex)
                        aLonOnEdge = aLongitudeVertex_180[aVertexOnEdgeIndex-1]
                        aLatOnEdge = aLatitudeVertex[aVertexOnEdgeIndex-1]
                        print("Edge ID:", aEdgeID[iEdge])
                        print("Vertices on edge:")
                        npoints = len(aLonOnEdge)
                        for iPoint in range(npoints):
                            x1 = aLonOnEdge[iPoint]
                            y1 = aLatOnEdge[iPoint]
                            print("Point:", x1, y1)
                        pass


                #vertices on cell
                aVertexOnCellIndex = aVertexOnCell[lCellIndex,:]
                pass
            else:
                #print("Polygon is valid...")
                pass



    return
#create a main call
if __name__ == '__main__':

    iFlag_plot_in = 1
    sFolder_out = '/qfs/people/liao313/workspace/python/hexwatershed_utility/figures/mpas'
    sPath_to_find = '/compyfs/liao313/04model/pyflowline/global/pyflowline20250101010/jigsaw/'
    #find under this path for each folder
    aFolder = os.listdir(sPath_to_find)
    for sFolder in aFolder:
        #get the length of the folder
        nLength = len(sFolder)
        if nLength > 3 and sFolder[0:3] == 'tmp':
            #this is the folder
            sFilename_mpas_mesh_netcdf_base = os.path.join(sPath_to_find, sFolder, 'mesh_out.nc')
            exam_mpas_base_mesh(sFilename_mpas_mesh_netcdf_base, iFlag_plot_in = None, sFolder_out = sFolder_out)





