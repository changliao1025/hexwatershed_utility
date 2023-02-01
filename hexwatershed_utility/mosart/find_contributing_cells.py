import json
import  sys
import getpass
from datetime import datetime

import numpy as np

from netCDF4 import Dataset
from pyearth.gis.gdal.read.gdal_read_envi_file import gdal_read_envi_file_multiple_band

from pyflowline.algorithms.auxiliary.gdal_functions import calculate_distance_based_on_lon_lat



def find_contributing_cells(aLongitude_in, aLatitude_in, aCellID, aCellID_downslope, dLongitude, dLatitude):
    """
    Should the data in 2D? or MPAS based?
    """

    #find the cell what matches with drainage area?
   
       
    aCellID_downslope=np.array(aCellID_downslope)
    nCell = aCellID_downslope.shape[0]

    aDistance = np.full(nCell, -9999, dtype=float)

    #calculate the distance betwen the cell with other cells
    #this method is only approximate
    for i in range(nCell):
        #this is a very simple function 
        dummy = np.power(aLongitude_in[i]- dLongitude, 2) + np.power(aLatitude_in[i] - dLatitude ,2)
        aDistance[i] = np.sqrt( dummy )
    
    aIndex_distance = np.argsort(aDistance)
    #is it possible to have multipl nearest?
    lCellIndex_nearest= aIndex_distance[0]
    lCellID_nearest = aCellID[ lCellIndex_nearest ]
    
    nSearch = 1 
    aCellID_contribution = [lCellID_nearest]  
    aCellIndex_contribution = [lCellIndex_nearest]    

    for iSearch in range(nSearch):
        iIndex_dummy = aIndex_distance[iSearch]
        lCellID = aCellID[iIndex_dummy]
          

        aCellID_downslope_table = [lCellID]
        aCellIndex_downslope= [iIndex_dummy]
        iFlag_finished = 0
        while (iFlag_finished != 1):
            aCellID_downslope_current= list()
            nDownslope = len(aCellID_downslope_table)
            for i in range(nDownslope ):
                lCellID_dummy = aCellID_downslope_table[i]
                dummy_index0 = np.where( aCellID_downslope ==  lCellID_dummy )
                dummy_index=dummy_index0[0]
                nUpslope = len(dummy_index) 
                if nUpslope> 0:
                    for j in range(nUpslope):
                        lCellID_dummy = aCellID[dummy_index[j]]
                        aCellIndex_contribution.append( dummy_index[j] )
                        aCellID_contribution.append( lCellID_dummy )
                        aCellID_downslope_current.append(lCellID_dummy )
            
            if len(aCellID_downslope_current) > 0:
                aCellID_downslope_table = aCellID_downslope_current 
            else:
                iFlag_finished = 1
            
            pass
                
      
    return lCellID_nearest,lCellIndex_nearest, aCellID_contribution, aCellIndex_contribution

