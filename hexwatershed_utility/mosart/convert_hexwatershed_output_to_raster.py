import numpy as np
from osgeo import gdal, osr, ogr
def convert_hexwatershed_output_to_raster(sMesh_type_in, 
sFilename_json_in, sFilename_spatial_reference_in ):

    #first we need to read the json file

    #this function should only support the latlon and regular grid
    #this will also affects the coordiante system 
    #for latlon, the coordinate system is WGS84
    #for square grid, the coordinate system should be defined by the user

    if sMesh_type_in == 'latlon':
        #define the wgs 
        pSpatialRef = osr.SpatialReference()
        pSpatialRef.ImportFromEPSG(4326)
        pass
    else:
        if sMesh_type_in == 'square':

            pass
        else:
            pass

  