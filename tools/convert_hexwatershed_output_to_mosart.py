import json
import os, sys
from netCDF4 import Dataset
import numpy as np
from datetime import datetime
import scipy
from scipy.io import netcdf
import getpass



def find_contributing_cells(aLongitude_in, aLatitude_in, aCellID, aCellID_downslope, dLongitude, dLatitude):
    """
    Should the data in 2D? or MPAS based?
    """

    #find the cell what matches with drainage area?
   
       
    aCellID_downslope=np.array(aCellID_downslope)
    nCell = aCellID_downslope.shape()
    aDistance = np.full(nCell, -9999, dtype=float)

    #calculate the distance betwen the cell with other cells
    #this method is only approximate
    for i in range(nCell):
        #this is a very simple function 
        dummy = np.power(aLongitude_in(i)- dLongitude, 2) + np.power(aLatitude_in[i] - dLatitude ,2)
        aDistance[i] = np.sqrt( dummy )
    
    aIndex_distance = np.argsort(aDistance)
    #is it possible to have multipl nearest?
    lCellIndex_nearest= aIndex_distance[0]
    lCellID_nearest = aCellID[ lCellIndex_nearest ]
    
    nSearch = 1    

    for iSearch in range(nSearch):
        iIndex_dummy = aIndex_distance[iSearch]
        lCellID = aCellID[iIndex_dummy]
        aCellID_contribution = list()    
        aCellIndex_contribution = list()   

        aCellID_downslope_table = [lCellID]
        aCellIndex_downslope= [iIndex_dummy]
        iFlag_finished = 0
        while (iFlag_finished != 1):
            aCellID_downslope_current= list()
            nDownslope = len(aCellID_downslope_table)
            for i in range(nDownslope ):
                dummy_index = np.where( aCellID_downslope ==  aCellID_downslope_table[i] )
                if len(dummy_index) > 0:
                    for j in dummy_index:
                        lCellID_dummy = aCellID[dummy_index[j]]
                        aCellIndex_contribution.append( dummy_index[j] )
                        aCellID_contribution.append( lCellID_dummy )
                        aCellID_downslope_current.append(lCellID_dummy )
            pass
            aCellID_downslope_table = aCellID_downslope_current 
                
      
    return lCellID_nearest,lCellIndex_nearest, aCellID_contribution, aCellIndex_contribution


def get_geometry(aLongitude_in, aLatitude_in, aCellID, aCellID_downslope, aArea, pWidth_in = None, pDepth_in = None):
    """
    This is the function to retrive the river width, depth and 
    """
    #generate the mesh 
    nCell = len(aLongitude_in)

    aLongitude0 = np.arange(-179.75,179.75,0.5)
    aLatitude0 = np.arange(-59.75,89.75, 0.5)
    aLongitude, aLatitude = np.meshgrid(aLongitude0, aLatitude0)

    aLongitude = np.transpose(aLongitude)
    aLatitude = np.transpose(aLatitude)

    if pWidth_in is None:
        pWidth = 7.2
    else:
        pWidth = pWidth_in
    if pDepth_in is None:
        pDepth = 0.27
    else:
        pDepth = pDepth_in     
    

    #initialize runoff and discharge
    aRunoff = np.full( (nCell, 31 * 365), -9999, dtype = float)
    aDischarge = np.full( (nCell, 31 * 365), -9999, dtype = float)
    AMF = np.full( (nCell, 31), -9999, dtype = float)
    k = 1
    print('Generating nearest neighbour mapping...')
    aIndexNearest = np.full((nCell) , -9999, dtype = int)
  
    for i in range( nCell):
        dummy0 = np.power( aLongitude - aLongitude_in(i), 2) 
        dummy1 = np.power( aLatitude - aLatitude_in(i), 2)

        aDistance = np.sqrt( dummy0 + dummy1 )
        
        dummy_index = np.where(aDistance == np.min(aDistance))

       
        if len(dummy_index) > 0:           
            index = dummy_index[0]       
            aIndexNearest[i] = index
    

    print('Reading daily aRunoff...\n')
    sWorkspace_runoff=''

    k=0
    for i in range (1979,2009,1):        
        si = "{:04d}".format(i)
        for j in range (1,365,1):
            sj = "{:04d}".format(j)
            sFilename = sWorkspace_runoff + '/' + 'RUNOFF05_' + si + '_' + sj + '.mat'
            #read mat uning 
            aData =  scipy.io.loadmat(sFilename)
            aRunoff[:, k] = aData[aIndexNearest]
            k = k + 1
       
    #find contributing cells
    dummy = np.full( nCell, -9999, dtype = 1 )
    aCellIndex_all= list()
    aCellID_all= list()
    aCellID_contribution_all= list()
    aCellIndex_contribution_all= list()
    print('Searching for contribuing area...\n')

    for i in range(0, nCell,1):      
        lCellID_nearest, lCellIndex_nearest, aCellID_contribution, aCellIndex_contribution = find_contributing_cells(aLongitude_in, aLatitude_in, aCellID, aCellID_downslope, \
            aArea,  aLongitude_in[i], aLatitude_in[i])   
        aCellIndex_all.append(lCellIndex_nearest)
        aCellID_all.append(lCellID_nearest)
        aCellIndex_contribution_all.append(aCellIndex_contribution)

    print('Mapping Runoff to Discharge...\n')

    #sun up runoff to get discrage
    for i in range(1, nCell,1):      
        dummy3 = aCellID_contribution_all[i]
        dummy0 = aRunoff[dummy3] 
        dummy1 = aArea[i]
        dummy2 = dummy0 * dummy1
        aDischarge[i, :] = np.sum( dummy2  / 1000 / (3 * 60 * 60))


    for i in range (1,31,1):
        tmp = aDischarge[:, (i - 1) * 365 + 1:i * 365]
        AMF[:, i] = max(tmp, [], 2)
   

    flood_2yr = np.percentile(AMF, 50, 2)

    rwid = aw * np.power(flood_2yr, 0.52)
    rdep = ad * np.power(flood_2yr, 0.31)

    return rwid, rdep, flood_2yr


def create_unstructure_domain_file_1d(aLon_region, aLat_region, \
   aLonV_region, aLatV_region,  
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
    var = dict()
    aVariable = ['area','frac','mask','aLongitude_in','xv','aLatitude_in','yv']
    nVariable = len(aVariable)
    aUnit = ['area','frac','mask','aLongitude_in','xv','aLatitude_in','yv']
    aLongName = ['area','frac','mask','aLongitude_in','xv','aLatitude_in','yv']
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
    #for varname in ncid_inq.variables:
    for i in range(nVariable):    
        varname = aVariable[i]
        if varname == 'aLongitude_in':
            data = aLon_region
        elif varname == 'aLatitude_in':
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
                
                pass
        
        var[varname][:] = data 


 
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
        lID =0 
        for i in range(ncell):
            pcell = data[i]
            aID.append(lID)
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
            aRdep.append(1.0)     
            aRwid.append(1.0)
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
    
    nCell = len(aCellID)
    aCellID = np.array(aCellID)
    aCellID_downslope=np.array(aCellID_downslope)
    #convert to numpy array
    aArea=np.array(aArea)
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
        
        
            index = np.reshape(dummy_index, 1)[0]
            aDnID.append(index + 1)
       
            
    aDnID = np.array(aDnID)

    
    aDomainfrac=np.array(aDomainfrac)
    aElevation=np.array(aElevation)
    
   
    aRlen=np.array(aFlowline_length)
    aLon  =np.array(aLongitude)
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
    aLonV_region, aLatV_region,     sFilename_mosart_domain_out)  

    return


if __name__ == '__main__':

    sFilename_json_in='/compyfs/liao313/04model/pyhexwatershed/sag/pyhexwatershed20220607001/hexwatershed/hexwatershed.json'
    sFilename_mpas_in='/people/liao313/workspace/python/pyhexwatershed_icom/data/sag/input/lnd_mesh.nc'
    sFilename_mosart_parameter_in = '/compyfs/inputdata/rof/mosart/MOSART_Global_half_20210616.nc'
    sFilename_mosart_parameter_out = 'mosart_sag_parameter.nc'
    sFilename_mosart_domain_out = 'mosart_sag_domain.nc'
    convert_hexwatershed_json_to_mosart_netcdf(sFilename_json_in, \
        sFilename_mpas_in, \
            sFilename_mosart_parameter_in,
            sFilename_mosart_parameter_out,\
            sFilename_mosart_domain_out)