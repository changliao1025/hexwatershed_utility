import os
import sys
import math
import numpy as np
import json
from osgeo import  osr, gdal
from pyearth.gis.gdal.read.raster.gdal_read_geotiff_file import gdal_read_geotiff_file
from pyearth.gis.gdal.write.raster.gdal_write_geotiff_file import gdal_write_geotiff_file
from pyearth.gis.spatialref.reproject_coodinates import reproject_coordinates, reproject_coordinates_batch

def index_to_row_col(index, num_columns):
    index -= 1  # Adjust for 1-based indexing
    row = index // num_columns + 1
    col = index % num_columns + 1
    return row, col

def calculate_flow_direction_index_from_row_column(lRow_index, lColumn_index, lRow_index1, lColumn_index1):
    row_diff = lRow_index1 - lRow_index
    col_diff = lColumn_index1 - lColumn_index

    if row_diff == 1 and col_diff == 0:
        return 7  # South
    elif row_diff == 1 and col_diff == -1:
        return 6  # SouthWest
    elif row_diff == 0 and col_diff == -1:
        return 5  # West
    elif row_diff == -1 and col_diff == -1:
        return 4  # NorthWest
    elif row_diff == -1 and col_diff == 0:
        return 3  # North
    elif row_diff == -1 and col_diff == 1:
        return 2  # NorthEast
    elif row_diff == 0 and col_diff == 1:
        return 1  # East
    elif row_diff == 1 and col_diff == 1:
        return 8  # SouthEast
    else:
        return -1  # No direction

def convert_hexwatershed_output_to_raster(sMesh_type_in,
                                          sFilename_json_in,
                                          sFilename_spatial_reference_in):

    # first we need to read the json file

    # this function should only support the latlon and regular grid
    # this will also affects the coordiante system
    # for latlon, the coordinate system is WGS84
    # for square grid, the coordinate system should be defined by the user
    pSpatialRef_source = osr.SpatialReference()
    pSpatialRef_source.ImportFromEPSG(4326)
    if sMesh_type_in == 'latlon':
        # define the wgs
        
        pass
    else:
        if sMesh_type_in == 'square':
            # what kind of project we should use? a raster based or a vector based?
            # in general, because the square case should be used from a DEM, then it is more reasonable to use a raster based projection
            # call the pyearth api to read the spatial reference

            dummy = gdal_read_geotiff_file(sFilename_spatial_reference_in)
            pSpatialRef_target = dummy['spatialReference']
            #nrow = dummy['nrow']
            #ncolumn = dummy['ncolumn']
            dPixelWidth = dummy['pixelWidth']
            dPixelHeight = dummy['pixelHeight']
            #dOriginX_in = dummy['originX']
            #dOriginY_in = dummy['originY']


            pass
        else:
            pass

    # read the json file
    #we will first generate DEM, slope and flow accumulation, and mask

    
    aLongitude = list()
    aLatitude = list()
    aCellID = dict()
    

    print('start calculating domain')
    sys.stdout.flush()
    with open(sFilename_json_in) as json_file:
        data = json.load(json_file)      
        ncell = len(data)
        lCellIndex = 0
        for i in range(ncell):
            pcell = data[i]            
            lCellID = pcell['lCellID']
            dLongitude = pcell['dLongitude_center_degree']
            dLatitude = pcell['dLatitude_center_degree']
            aLongitude.append(dLongitude)
            aLatitude.append(dLatitude)
            #aCellID.append(lCellID)
            aCellID[lCellID] = lCellIndex
            lCellIndex = lCellIndex + 1
    
    #now reproject them back to the original coordinate system
    aX, aY = reproject_coordinates_batch(aLongitude, aLatitude, pSpatialRef_source, pSpatialRef_target)

    #get min and max of x and y
    dX_min = min(aX)
    dX_max = max(aX)
    dY_min = min(aY)
    dY_max = max(aY)
    dPixelHeight1 = np.abs(dPixelHeight)

    nrow = int((dY_max - dY_min)/dPixelHeight1) + 1
    ncolumn = int((dX_max - dX_min)/dPixelWidth) + 1

    dOriginX_in = dX_min - dPixelWidth/2
    dOriginY_in = dY_max + dPixelHeight1/2

    aDEM = np.full((nrow, ncolumn), -9999.0, dtype=np.float32)
    aSlope = np.full((nrow, ncolumn), -9999.0, dtype=np.float32)
    aFlow_accumulation = np.full((nrow, ncolumn), -9999, dtype=int)
    aMask = np.full((nrow, ncolumn), -9999, dtype=int)
    aFlow_direction = np.full((nrow, ncolumn), 0, dtype=int)       

    print('start calculating variables')
    sys.stdout.flush()
    with open(sFilename_json_in) as json_file:
        data = json.load(json_file)      
        ncell = len(data)
       
        for i in range(ncell):
            pcell = data[i]            
           
            #we now use the longitude and latitude to determine the row and column index

            dLongitude = pcell['dLongitude_center_degree']
            dLatitude = pcell['dLatitude_center_degree']

            dX, dY = reproject_coordinates(dLongitude, dLatitude, pSpatialRef_source, pSpatialRef_target)    

            lRow_index = int(np.round((dY_max - dY)/dPixelHeight1) )
            lColumn_index = int(np.round((dX - dX_min)/dPixelWidth)  )
            #print(lRow_index, lColumn_index)  
                
            aDEM[lRow_index, lColumn_index] = float(pcell['dElevation'])
            aSlope[lRow_index, lColumn_index] = float(pcell['dSlope_between'])
            aFlow_accumulation[lRow_index, lColumn_index] = int(pcell['dDrainage_area'] / (dPixelHeight1 * dPixelWidth))

            if int(pcell['lStream_segment']) >0 :
                aMask[lRow_index, lColumn_index] = 1 #stream
            else:
                aMask[lRow_index, lColumn_index] = 0 #non-stream

            #now calcualte flow direction 
            lCellID_downslope = pcell['lCellID_downslope']
            if lCellID_downslope != -1:
                lCellIndex = aCellID[lCellID_downslope]
                pcell1 = data[lCellIndex]
                dLongitude1 = pcell1['dLongitude_center_degree']
                dLatitude1 = pcell1['dLatitude_center_degree']
                dX1, dY1 = reproject_coordinates(dLongitude1, dLatitude1, pSpatialRef_source, pSpatialRef_target)
                lRow_index1 = int( np.round((dY_max - dY1)/dPixelHeight1))
                lColumn_index1 = int(np.round((dX1 - dX_min)/dPixelWidth) )   
                #use row and columne index to determine the flow direction
                iFlow_direction = calculate_flow_direction_index_from_row_column(lRow_index, lColumn_index, lRow_index1, lColumn_index1)
                aFlow_direction[lRow_index, lColumn_index] = iFlow_direction   
            else:
                aFlow_direction[lRow_index, lColumn_index] = 5     
    

    #then we will generate the flow direction, which requires more complex algorithm

    #can we use the dem orginal? yes, because there are the same 

    print('start saving file')
    #flush io
    sys.stdout.flush()
    #write each individual raster file to the disk
    sWorkspace_out = '/compyfs/liao313/04model/pyhexwatershed/k34/pyhexwatershed20231001003/hexwatershed/00000001'

    #define output file name 
    sFilename_dem = os.path.join(sWorkspace_out, 'dem.tif')
    
    sFilename_slope = os.path.join(sWorkspace_out, 'slope.tif')
    sFilename_flow_accumulation = os.path.join(sWorkspace_out, 'flow_accumulation.tif')
    sFilename_mask = os.path.join(sWorkspace_out, 'mask.tif')
    sFilename_flow_direction = os.path.join(sWorkspace_out, 'flow_direction.tif')

    dMissing_value_in = -9999.0
    gdal_write_geotiff_file(sFilename_dem,
                            aDEM,
                            dPixelWidth,
                            dPixelHeight,
                            dOriginX_in,
                            dOriginY_in,
                            dMissing_value_in,
                            pSpatialRef_target)
    
    dMissing_value_in = -9999.0
    gdal_write_geotiff_file(sFilename_slope,
                            aSlope,
                            dPixelWidth,
                            dPixelHeight,
                            dOriginX_in,
                            dOriginY_in,
                            dMissing_value_in,
                            pSpatialRef_target)
    
    dMissing_value_in = -9999
    gdal_write_geotiff_file(sFilename_flow_accumulation,
                            aFlow_accumulation,
                            dPixelWidth,
                            dPixelHeight,
                            dOriginX_in,
                            dOriginY_in,
                            dMissing_value_in,
                            pSpatialRef_target,
                            datatype=gdal.GDT_Int32 )
    
    dMissing_value_in = -9999
    gdal_write_geotiff_file(sFilename_mask,
                            aMask,
                            dPixelWidth,
                            dPixelHeight,
                            dOriginX_in,
                            dOriginY_in,
                            dMissing_value_in,
                            pSpatialRef_target,
                            datatype=gdal.GDT_Int32   )
    

    dMissing_value_in = 0
    gdal_write_geotiff_file(sFilename_flow_direction,
                            aFlow_direction,
                            dPixelWidth,
                            dPixelHeight,
                            dOriginX_in,
                            dOriginY_in,
                            dMissing_value_in,
                            pSpatialRef_target,
                            datatype=gdal.GDT_Int32 )


