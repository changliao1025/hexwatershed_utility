
import  sys


import numpy as np

from netCDF4 import Dataset
from pyearth.gis.gdal.read.gdal_read_envi_file import gdal_read_envi_file_multiple_band
from hexwatershed_utility.mosart.find_contributing_cells import find_contributing_cells
def get_geometry(aLongitude_in, aLatitude_in, aCellID, aCellID_downslope, aArea, pWidth_in = None, pDepth_in = None):
    """
    This is the function to retrive the river width, depth and 

    Args:
        aLongitude_in (_type_): _description_
        aLatitude_in (_type_): _description_
        aCellID (_type_): _description_
        aCellID_downslope (_type_): _description_
        aArea (_type_): _description_
        pWidth_in (_type_, optional): _description_. Defaults to None.
        pDepth_in (_type_, optional): _description_. Defaults to None.

    Returns:
        _type_: _description_
    """
    #generate the mesh 
    nCell = len(aLongitude_in)
    iFlag_debug = 0
    if iFlag_debug==1:
        nCell=100

    dResolution_runoff = 0.5
    aLongitude0 = np.arange(-179.75,180,dResolution_runoff)
    aLatitude0 = np.arange(-59.75,90, dResolution_runoff)
    aLongitude, aLatitude = np.meshgrid(aLongitude0, aLatitude0)
    nrow, ncolumn = aLongitude.shape


    if pWidth_in is None:
        pWidth = 7.2
    else:
        pWidth = pWidth_in
    if pDepth_in is None:
        pDepth = 0.27
    else:
        pDepth = pDepth_in     
    

    #initialize runoff and discharge
    nyear = 2009-1978
    nday = nyear * 365
    aRunoff = np.full( (nday, nCell ), 0.0, dtype = float)
    aDischarge = np.full( ( nday, nCell), 0.0, dtype = float)
    AMF = np.full( (nyear, nCell ), 0.0, dtype = float)
    k = 1
    print('Generating nearest neighbour mapping...')
    #aIndexNearest = np.full((nCell) , -9999, dtype = int)
    
    iFlag_read_runoff = 1
    if iFlag_read_runoff ==1:
        aIndexNearest=list()
        for i in range(nCell):        
            dummy0 = np.power( aLongitude - aLongitude_in[i], 2) 
            dummy1 = np.power( aLatitude - aLatitude_in[i], 2)
            aDistance = np.sqrt( dummy0 + dummy1 )
            #aDistance = [ calculate_distance_based_on_lon_lat(aLongitude_in[i], aLatitude_in[i], aLongitude[j,k],aLatitude[j,k]) \
            #    for j in range(nrow) for k in range(ncolumn) ]
            distance_min = np.min(aDistance)        
            dummy_index = np.where(aDistance == distance_min)
            row_index = dummy_index[0]
            column_index=dummy_index[1]
            if len(row_index) ==1:     
                aIndexNearest.append([row_index[0],column_index[0]])
            else:
                #more than one
                aIndexNearest.append([row_index[0],column_index[0]])
                pass
            
        aIndexNearest=np.array(aIndexNearest)
        aIndexNearest=np.transpose(aIndexNearest)
        aIndexNearest=tuple(map(tuple, aIndexNearest))
        
        print('Reading daily Runoff...\n')
        sWorkspace_runoff='/compyfs/liao313/00raw'

        #k=0
        #result_ids = []
        sFilename = sWorkspace_runoff + '/' + 'runoff.dat'
        dummy0 = gdal_read_envi_file_multiple_band(sFilename)
        dummy = dummy0[0]
        for i in range (nday): 
            d = dummy[ i,:]
            c = d[aIndexNearest]
            c[np.where(c==-9999.0)]=0.0
            aRunoff[i,:] = c
        
       
    #find contributing cells
    
    aCellIndex_all= list()
    aCellID_all= list()
    aCellID_contribution_all= list()
    aCellIndex_contribution_all= list()

    iFlag_search_contribution_cell=1
    if iFlag_search_contribution_cell ==1:
        print('Searching for contribuing area...\n')
        sys.stdout.flush()
        for i in range(0, nCell,1):    
            lCellID = aCellID[i]          
            aCellID_contribution, aCellIndex_contribution = find_contributing_cells(aCellID, aCellID_downslope, lCellID)

            aCellID_contribution_all.append(aCellID_contribution)
            #aCellID_all.append(lCellID_nearest)
            aCellIndex_contribution_all.append(aCellIndex_contribution)
    
    #aCellIndex_all = ray.get(aCellIndex_all)
    #aCellID_all = ray.get(aCellID_all) 
    #aCellIndex_contribution_all = ray.get(aCellIndex_contribution_all) 

    print('Mapping Runoff to Discharge...\n')
    sys.stdout.flush()

    #sun up runoff to get discrage
    for i in range(nCell):      

        lCellID = aCellID[i]        
        if lCellID ==152663:
            print('debug')    
        if lCellID ==150992:
            print('debug') 
        dummy3 = aCellIndex_contribution_all[i]
        dummy4 = np.array(dummy3)
        dummy0 = aRunoff[:,dummy4] 
        dummy1 = aArea[dummy4]
        dummy2 = dummy0 * dummy1
        dummy4 = np.sum( dummy2, 1 )
        #unit conersion
        aDischarge[:, i] = dummy4  / 1000 / (3 * 60 * 60)


    for i in range(0,nyear,1):
        dummy_range = np.arange(i * 365, (i+1) * 365 , 1)
        tmp = aDischarge[dummy_range,: ]
        AMF[i,:] = np.max(tmp, 0)
   

    aFlood_2yr_out = np.percentile(AMF, 50, 0)

    aWidth_out = pWidth * np.power(aFlood_2yr_out, 0.52)
    aDepth_out = pDepth * np.power(aFlood_2yr_out, 0.31)

    return aWidth_out, aDepth_out, aFlood_2yr_out


