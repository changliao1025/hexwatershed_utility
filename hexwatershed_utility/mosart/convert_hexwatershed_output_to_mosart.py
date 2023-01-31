import json
import os, sys
from netCDF4 import Dataset
import numpy as np
from datetime import datetime

from scipy.io import netcdf
import getpass

from pyflowline.algorithms.auxiliary.gdal_functions import calculate_distance_based_on_lon_lat
from pyearth.gis.gdal.read.gdal_read_envi_file import gdal_read_envi_file_multiple_band



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


def get_geometry(aLongitude_in, aLatitude_in, aCellID, aCellID_downslope, aArea, pWidth_in = None, pDepth_in = None):
    """
    This is the function to retrive the river width, depth and 
    """
    #generate the mesh 
    nCell = len(aLongitude_in)
    iFlag_debug = 0
    if iFlag_debug==1:
        nCell=100

    aLongitude0 = np.arange(-179.75,180,0.5)
    aLatitude0 = np.arange(-59.75,90, 0.5)
    aLongitude, aLatitude = np.meshgrid(aLongitude0, aLatitude0)
    nrow, ncolumn = aLongitude.shape

    #aLongitude = np.transpose(aLongitude)
    #aLatitude = np.transpose(aLatitude)

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
    #aCellID_contribution_all= list()
    aCellIndex_contribution_all= list()

    iFlag_search_contribution_cell=1
    if iFlag_search_contribution_cell ==1:
        print('Searching for contribuing area...\n')
        sys.stdout.flush()
        for i in range(0, nCell,1):      
            lCellID_nearest, lCellIndex_nearest, aCellID_contribution, aCellIndex_contribution = find_contributing_cells(aLongitude_in, aLatitude_in, aCellID, aCellID_downslope, \
              aLongitude_in[i], aLatitude_in[i])   

            aCellIndex_all.append(lCellIndex_nearest)
            aCellID_all.append(lCellID_nearest)
            aCellIndex_contribution_all.append(aCellIndex_contribution)
    
    #aCellIndex_all = ray.get(aCellIndex_all)
    #aCellID_all = ray.get(aCellID_all) 
    #aCellIndex_contribution_all = ray.get(aCellIndex_contribution_all) 

    print('Mapping Runoff to Discharge...\n')
    sys.stdout.flush()

    #sun up runoff to get discrage
    for i in range(0, nCell,1):      
        dummy3 = aCellIndex_contribution_all[i]
        dummy4 = np.array(dummy3)
        dummy0 = aRunoff[:,dummy4] 
        dummy1 = aArea[i]
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


def create_unstructure_domain_file_1d(aLon_region, aLat_region, \
   aLonV_region, aLatV_region,  aArea,
    sFilename_domain_file_out):

    #sFilename_domain_file_out = '%s/domain_%s_%s.nc' % \
    #            (out_netcdf_dir, clm_usrdat_name, datetime.now().strftime('c%-y%m%d'))
    print('  domain: ' + sFilename_domain_file_out)

    # Check if the file is available   
    
    pDatasets_out = Dataset(sFilename_domain_file_out, 'w',format="NETCDF3_CLASSIC")

    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #
    #                           Define dimensions
    #
    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
    ni, nv = aLonV_region.shape
    nj = 1

  
    dimname = 'ni'
    pDatasets_out.createDimension(dimname,ni)
    dimname = 'nj'
    pDatasets_out.createDimension(dimname,nj)
    dimname = 'nv'
    pDatasets_out.createDimension(dimname,nv)

    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #
    #                           Define variables
    #
    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    aDimension_list1=list()
    aDimension_list1.append('nj')
    aDimension_list1.append('ni')
    aDimension_tuple1 = tuple(aDimension_list1)
    aDimension_list2=list()
    aDimension_list2.append('nj')
    aDimension_list2.append('ni')
    aDimension_list2.append('nv')
    aDimension_tuple2 = tuple(aDimension_list2)
    
    aVariable = ['area','frac','mask','xc','xv','yc','yv']
    nVariable = len(aVariable)
    aUnit = ['area','frac','mask','xc','xv','yc','yv']
    aLongName = ['area','frac','mask','xc','xv','yc','yv']
    for i in range(nVariable):
        varname = aVariable[i]
        if varname == 'xv' or varname == 'yv':
            dtype = float
            dims  = aDimension_tuple2
        else:
            dtype = float
            dims  = aDimension_tuple1

        pVar = pDatasets_out.createVariable(varname, dtype, dims, fill_value=-9999)
        
        pVar.setncatts( { '_FillValue': -9999 } )
        pVar.setncatts( { 'unit':  aUnit[i] } )
        pVar.setncatts( { 'long name': aLongName[i] } )    

    user_name = getpass.getuser()
    setattr(pDatasets_out,'Created_by',user_name)
    setattr(pDatasets_out,'Created_on',datetime.now().strftime('%c'))

    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #
    #                           Copy variables
    #
    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    for sKey, aValue in pDatasets_out.variables.items():
    #for i in range(nVariable):    
        varname = sKey#aVariable[i]
        if varname == 'xc':
            data = aLon_region
        elif varname == 'yc':
            data = aLat_region
        elif varname == 'xv':
            data = aLonV_region
        elif varname == 'yv':
            data = aLatV_region
        elif varname == 'mask':
            data = np.ones( (1,len(aLon_region)) )
        elif varname == 'frac':
            data = np.ones( (1,len(aLon_region)) )
        elif varname == 'area':
            if aLonV_region.shape[1] == 3:
                ax = aLonV_region[:,0]
                ay = aLatV_region[:,0]
                bx = aLonV_region[:,1]
                by = aLatV_region[:,1]
                cx = aLonV_region[:,2]
                cy = aLatV_region[:,2]

                data = 0.5*(ax*(by-cy) + bx*(cy-ay) + cx*(ay-by))
            elif aLonV_region.shape[1] == 4:
                data = (aLonV_region[:,0] - aLonV_region[:,1]) * (aLatV_region[:,0] - aLatV_region[:,2])
            else:
                #raise NameError('Added area computation')
                #use different method to get area
                radius= 6378137.0                      
                dummy_data = np.array(aArea )
                data  = dummy_data / ( 4*np.pi*(radius**2) )

                pass
        
        aValue[:] = data 


 
    pDatasets_out.close()

    return sFilename_domain_file_out


def convert_hexwatershed_json_to_mosart_netcdf(sFilename_json_in, \
    sFilename_mpas_in, \
    sFilename_mosart_parameter_in, \
        sFilename_mosart_parameter_out, \
            sFilename_mosart_domain_out):
    #open json and read data
    aID =list()
    aDnID=list()
    aCellID = list()
    aCellID_downslope=list()
    aElevation=list()
    aDrainage=list()
    aSlope=list()
    aFlowline_length=list()
    aLongitude  =list()
    aLatitude =list()
    aArea=list()
    aDomainfrac=list()
    aFdir =list()
    aGxr =list()

    aHslp =list()
    aID=list()
    aLat=list()
    aLatixy=list()
    aLon=list()
    aLongxy=list()
    aNh=list()
    aNr=list()
    aNt=list()
    aRdep=list()
    aRlen=list()
    aRslp=list()
    aRwid=list()
    aTslp=list()
    aTwid=list()
    aLatV_region=list()
    aLonV_region=list()
    with open(sFilename_json_in) as json_file:
        data = json.load(json_file)      
        ncell = len(data)
        lID = 1 
        for i in range(ncell):
            pcell = data[i]
            aID.append(lID)
            lID = lID + 1
            dArea = float(pcell['Area'])
            aArea.append(dArea)
            aLongitude.append(float(pcell['dLongitude_center_degree']))
            aLatitude.append(float(pcell['dLatitude_center_degree']))
            aCellID.append(int(pcell['lCellID']))
            aCellID_downslope.append(int(pcell['lCellID_downslope']))
            aElevation.append(float(pcell['Elevation']))
            dDrainage = float(pcell['DrainageArea'])
            aDrainage.append(dDrainage)
            aSlope.append(float(pcell['dSlope_between']))
            dFlowline_length = float(pcell['dLength_flowline'])
            aFlowline_length.append(dFlowline_length)
            aDomainfrac.append(1.0)
            aFdir.append(1)

            dGxr = dArea / dFlowline_length
            aGxr.append(dGxr)
               
            aNh.append(1.0)
            aNr.append(1.0)
            aNt.append(1.0)
            
            aTwid.append(1.0)

            
            dummy_vertex = pcell['vVertex']
            aVertex_lon = np.full(8, -9999, float)
            aVertex_lat = np.full(8, -9999, float)
            nVertex=len(dummy_vertex)
            for j in range(nVertex):
                aVertex_lon[j] = dummy_vertex[j]['dLongitude_degree']
                aVertex_lat[j] = dummy_vertex[j]['dLatitude_degree']

            aLonV_region.append(aVertex_lon)
            aLatV_region.append(aVertex_lat)
            
        pass
    
    aLongitude_in = np.array(aLongitude)
    aLatitude_in = np.array(aLatitude)
    nCell = len(aCellID)
    aCellID = np.array(aCellID)
    aCellID_downslope=np.array(aCellID_downslope)
    aArea=np.array(aArea)
    aRwid, aRdep, aFlood_2yr_out = get_geometry(aLongitude_in, aLatitude_in, aCellID, aCellID_downslope, aArea, pWidth_in = None, pDepth_in = None)
    
    
    
    #convert to numpy array
    
    aAreaTotal2=np.array(aDrainage)
    aAreaTotal= aAreaTotal2
    
    aID = np.array(aID)

    aLonV_region= np.array(aLonV_region)
    aLatV_region= np.array(aLatV_region)

    for i in range(nCell):
        lCellID = aCellID[i]
        lCellID_downslope = aCellID_downslope[i]
        if lCellID_downslope == -1:
            aDnID.append(-9999)
        else:
            dummy_index  = np.where(aCellID == lCellID_downslope)
            #index = np.reshape(dummy_index, 1)[0]
            dnID = aID[dummy_index[0]]
            aDnID.append(dnID[0])
       
    aDnID = np.array(aDnID)
    aDomainfrac=np.array(aDomainfrac)
    aElevation=np.array(aElevation)
    aRlen=np.array(aFlowline_length)
    aLon =np.array(aLongitude)
    aLat =np.array(aLatitude)
    aLongxy = aLon
    aLatixy = aLat


    aDomainfrac = np.array(aDomainfrac)
    aFdir=np.array(aFdir)
    aGxr=np.array(aGxr)
    aHslp=np.array(aSlope)     
    aNh=np.array(aNh)
    aNr=np.array(aNr)
    aNt=np.array(aNt)
    aRdep=np.array(aRdep)     
    aRwid=np.array(aRwid)
    aTwid=np.array(aTwid)  
    aRslp = aHslp
    aTslp = aHslp

    #open mpas netcdf 
    pDatasets_in = Dataset(sFilename_mosart_parameter_in)
    netcdf_format = pDatasets_in.file_format
    #output file
    pDatasets_out = Dataset(sFilename_mosart_parameter_out, "w", format=netcdf_format)
    aVariable = ['area','areaTotal','areaTotal2','dnID', 'domainfrac',\
        'fdir','frac','gxr','hslp','ID',\
            'lat','latixy','lon', 'longxy','nh', 'nr','nt',
            'rdep','rlen','rslp','rwid', 'tslp', 'twid']

    aUnit = ['','','','','','','','','','','','','','','','','','','','','','','','','','','','','','']
    aLongName = ['','','','','','','','','','','','','','','','','','','','','','','','','','','','','','']

    pDatasets_out.createDimension('gridcell', nCell)

    #write data out

    aData_out=list()
    aData_out.append(aArea)
    aData_out.append(aAreaTotal)
    aData_out.append(aAreaTotal2)
    aData_out.append(aDnID)
    aData_out.append(aDomainfrac)
    aData_out.append(aFdir)
    aData_out.append(aFdir)
    aData_out.append(aGxr)
    aData_out.append(aHslp)
    aData_out.append(aID)
    aData_out.append(aLat)
    aData_out.append(aLatixy)
    aData_out.append(aLon)
    aData_out.append(aLongxy)
    aData_out.append(aNh)
    aData_out.append(aNr)
    aData_out.append(aNt)
    aData_out.append(aRdep)
    aData_out.append(aRlen)
    aData_out.append(aRslp)
    aData_out.append(aRwid)
    aData_out.append(aTslp)
    aData_out.append(aTwid)

    #now export
    nVariable = len(aVariable)
    aDimension_list=list()

    aDimension_list.append('gridcell')
    aDimension_tuple = tuple(aDimension_list)
    for i in range(nVariable):   
        pVar = pDatasets_out.createVariable(aVariable[i], float, aDimension_tuple, fill_value=-9999  )
        aData_variable = aData_out[i]
        print(aVariable[i])
        pVar[:] = aData_variable
        pVar.setncatts( { '_FillValue': -9999 } )
        pVar.setncatts( { 'unit':  aUnit[i] } )
        pVar.setncatts( { 'long name': aLongName[i] } )  

    
    #create domain file
    
    aLon_region  = aLon
    aLat_region = aLat
  

  
    create_unstructure_domain_file_1d(aLon_region, aLat_region, \
    aLonV_region, aLatV_region,  aArea,   sFilename_mosart_domain_out)  

    return


if __name__ == '__main__':
    sRegion ='columbia'

    if sRegion == 'sag':

        sFilename_json_in='/compyfs/liao313/04model/pyhexwatershed/sag/pyhexwatershed20220607001/hexwatershed/hexwatershed.json'
        sFilename_mpas_in='/people/liao313/workspace/python/pyhexwatershed_icom/data/sag/input/lnd_mesh.nc'
        sFilename_mosart_parameter_in = '/compyfs/inputdata/rof/mosart/MOSART_Global_half_20210616.nc'
        sFilename_mosart_parameter_out = 'mosart_sag_parameter.nc'
        sFilename_mosart_domain_out = 'mosart_sag_domain.nc'
    else:
        if sRegion == 'susquehanna':
            sFilename_json_in='/compyfs/liao313/04model/pyhexwatershed/susquehanna/pyhexwatershed20221115001/hexwatershed/hexwatershed.json'
            sFilename_mpas_in='/qfs/people/liao313/workspace/python/pyhexwatershed_icom/data/susquehanna/input/lnd_cull_mesh.nc'
            sFilename_mosart_parameter_in = '/compyfs/inputdata/rof/mosart/MOSART_Global_half_20210616.nc'
            sFilename_mosart_parameter_out = 'mosart_susquehanna_parameter.nc'
            sFilename_mosart_domain_out = 'mosart_susquehanna_domain.nc'
        else:
            sFilename_json_in='/compyfs/liao313/04model/pyhexwatershed/columbia/pyhexwatershed20221115003/hexwatershed/hexwatershed.json'
            sFilename_mpas_in='/compyfs/liao313/00raw/mesh/global/lnd_mesh.nc'
            sFilename_mosart_parameter_in = '/compyfs/inputdata/rof/mosart/MOSART_Global_half_20210616.nc'
            sFilename_mosart_parameter_out = 'mosart_columbia_parameter.nc'
            sFilename_mosart_domain_out = 'mosart_columbia_domain.nc'
    
    convert_hexwatershed_json_to_mosart_netcdf(sFilename_json_in, \
        sFilename_mpas_in, \
            sFilename_mosart_parameter_in,
            sFilename_mosart_parameter_out,\
            sFilename_mosart_domain_out)