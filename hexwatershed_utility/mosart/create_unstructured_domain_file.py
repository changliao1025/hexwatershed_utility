
import getpass
from datetime import datetime
import numpy as np
from netCDF4 import Dataset

def create_unstructured_domain_file(aLon_region, aLat_region, \
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
    nv, ni = aLonV_region.shape
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

                #data = 0.5*(ax*(by-cy) + bx*(cy-ay) + cx*(ay-by))
            elif aLonV_region.shape[1] == 4:
                #data = (aLonV_region[:,0] - aLonV_region[:,1]) * (aLatV_region[:,0] - aLatV_region[:,2])
                
                pass
            else:
                #raise NameError('Added area computation')
                #use different method to get area
                radius= 6378137.0                      
                dummy_data = np.array(aArea ) #m^2
                data  = dummy_data / ( 4*np.pi*(radius**2) )

                pass
        
        aValue[:] = data 

    pDatasets_out.close()

    return sFilename_domain_file_out
