
from netCDF4 import Dataset #read netcdf
from osgeo import gdal, osr #the default operator
import numpy as np
from scipy.io import loadmat
from pyearth.system.define_global_variables import *    
from pyearth.gis.gdal.write.gdal_write_envi_file import gdal_write_envi_file_multiple_band
from pyearth.gis.gdal.write.gdal_write_geotiff_file import gdal_write_geotiff_file_multiple_band
sWorkspace_runoff='/qfs/people/xudo627/ming/runoff'
sWorkspace_runoff_out='/compyfs/liao313/00raw/'

    

def compress_runoff_data():
    ncolumn = 720
    nrow = 300
    dLon_min = -179.75
    dLat_max = 89.75
    dResolution_x =0.5
    dResolution_y =0.5
    pHeaderParameters = {}    
    pHeaderParameters['ncolumn'] = "{:0d}".format(ncolumn)
    pHeaderParameters['nrow'] = "{:0d}".format(nrow)
    pHeaderParameters['ULlon'] = "{:0f}".format(dLon_min-0.5 * dResolution_x)
    pHeaderParameters['ULlat'] = "{:0f}".format(dLat_max+0.5 * dResolution_y)
    pHeaderParameters['pixelSize'] = "{:0f}".format(dResolution_x)
    pHeaderParameters['nband'] = '1'
    pHeaderParameters['offset'] = '0'
    pHeaderParameters['data_type'] = '4'
    pHeaderParameters['bsq'] = 'bsq'
    pHeaderParameters['byte_order'] = '0'
    pHeaderParameters['missing_value'] = '-9999'
    iYear_end = 2009
    iYear_start = 1979
    nday = (iYear_end - iYear_start +1) * 365
    aGrid_stack= np.full((nday, nrow, ncolumn), -9999.0, dtype= float)
    k=0
    for i in range (1979,2009+1,1):        
        si = "{:04d}".format(i)
        for j in range (1,366,1):
            sj = "{:0d}".format(j)
            
            sFilename = sWorkspace_runoff + '/' + 'RUNOFF05_' + si + '_' + sj + '.mat'
            print(sFilename)
            dummy_data =  loadmat(sFilename)
            dummy = np.transpose(dummy_data['ro05'])
            nan_index = np.where(np.isfinite(dummy)==False)
            dummy[nan_index] = -9999
            aGrid_stack[k, :,:] = dummy
            k= k+1


    pSpatial = osr.SpatialReference()
    pSpatial.ImportFromEPSG(4326)
    sFilename_envi = sWorkspace_runoff_out + slash + 'runoff'  + sExtension_envi

    gdal_write_envi_file_multiple_band(sFilename_envi, aGrid_stack,\
        float(pHeaderParameters['pixelSize']),\
         float(pHeaderParameters['ULlon']),\
              float(pHeaderParameters['ULlat']),\
                  -9999.0, pSpatial)

if __name__ == '__main__':
    compress_runoff_data()