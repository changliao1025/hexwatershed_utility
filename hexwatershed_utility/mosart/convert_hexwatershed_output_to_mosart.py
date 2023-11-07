import os
from pathlib import Path
import json
import numpy as np
import netCDF4 as nc
from pye3sm.mesh.unstructured.e3sm_create_unstructured_domain_file_full import e3sm_create_unstructured_domain_file_full
from hexwatershed_utility.mosart.get_geometry import get_geometry
def convert_hexwatershed_json_to_mosart_netcdf(sFilename_json_in, \
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
    aRwid0=list()
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
            #qa for slope
            dSlope = float(pcell['dSlope_between'])
            if dSlope < 0.0001:
                dSlope = 0.0001
                pass
            
            aSlope.append(dSlope)
            dFlowline_length = float(pcell['dLength_flowline'])
            aFlowline_length.append(dFlowline_length)
            aDomainfrac.append(1.0)
            aFdir.append(1)

            #dGxr = dArea / dFlowline_length
            dGxr =  dFlowline_length / dArea
            aGxr.append(dGxr)
               
            aNh.append(0.1)
            aNr.append(0.05)
            aNt.append(0.05)            
            aTwid.append(10.0)
            
            dummy_vertex = pcell['vVertex']
            aVertex_lon = np.full(9, -9999, float)
            aVertex_lat = np.full(9, -9999, float)
            nVertex=len(dummy_vertex)
            for j in range(nVertex):
                aVertex_lon[j] = dummy_vertex[j]['dLongitude_degree']
                aVertex_lat[j] = dummy_vertex[j]['dLatitude_degree']

            aLonV_region.append(aVertex_lon)
            aLatV_region.append(aVertex_lat)
            
        pass
    
    aLongitude_in = np.array(aLongitude).reshape(ncell)
    aLatitude_in = np.array(aLatitude).reshape(ncell)
    nCell = len(aCellID)
    aCellID = np.array(aCellID).reshape(ncell)
    aCellID_downslope=np.array(aCellID_downslope).reshape(ncell)
    aArea=np.array(aArea).reshape(ncell)
    
    aRwid, aRdep, aFlood_2yr_out = get_geometry(aLongitude_in, aLatitude_in, aCellID, aCellID_downslope, aArea, pWidth_in = None, pDepth_in = None)
    
    #convert to numpy array
    
    aAreaTotal2=np.array(aDrainage)
    aAreaTotal= aAreaTotal2
    
    aID = np.array(aID)

    aLonV_region= np.array(aLonV_region) #.T
    aLatV_region= np.array(aLatV_region) #.T



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
    pDatasets_in = nc.Dataset(sFilename_mosart_parameter_in)
    netcdf_format = pDatasets_in.file_format
    #output file
    pDatasets_out = nc.Dataset(sFilename_mosart_parameter_out, "w", format=netcdf_format)
    aVariable = ['area','areaTotal','areaTotal2','dnID', 'domainfrac',\
        'fdir','frac','gxr','hslp','CellID', 'ID',\
            'lat','latixy','lon', 'longxy','nh', 'nr','nt',
            'rdep','rlen','rslp','rwid','rwid0', 'tslp', 'twid']

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
    aData_out.append(aCellID) #global cell ID
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
    aRwid0=aRwid * 5.0
    aData_out.append(aRwid0)
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

    aShape = aLon.shape
    if len(aShape)==1:
        nrow=aShape[0]
        ncolumn = 1
    else:
        nrow=aShape[0]
        ncolumn = aShape[1]

    
    
    aLon.shape=(nrow, 1)
    aLat.shape=(nrow, 1)
    
    nrow, nvertex = aLonV_region.shape
    aLonV_region.shape = (nrow, 1, nvertex)
    aLatV_region.shape = (nrow, 1, nvertex)

    

    e3sm_create_unstructured_domain_file_full(aLon, aLat,  aLonV_region, aLatV_region,  
                              sFilename_mosart_domain_out, aArea_in=aArea )
    #two options: simple ; non-simple



    #e3sm_create_unstructured_domain_file_simple(aLon, aLat,  aLonV_region, aLatV_region,  
    #                          sFilename_mosart_domain_out, aArea_in=aArea )

    sFolder = os.path.dirname(sFilename_mosart_domain_out)
    sBasename = Path(sFilename_mosart_domain_out).stem
    sFilname_new = sBasename +  '_simple.nc'
    sFilename_mosart_domain_out= os.path.join(sFolder , sFilename_mosart_domain_out) 
    


    
    return


if __name__ == '__main__':
    sRegion ='columbia'

    if sRegion == 'sag':

        sFilename_json_in='/compyfs/liao313/04model/pyhexwatershed/sag/pyhexwatershed20220607001/hexwatershed/hexwatershed.json'
        #sFilename_mpas_in='/people/liao313/workspace/python/pyhexwatershed_icom/data/sag/input/lnd_mesh.nc'
        sFilename_mosart_parameter_in = '/compyfs/inputdata/rof/mosart/MOSART_Global_half_20210616.nc'
        sFilename_mosart_parameter_out = 'mosart_sag_parameter.nc'
        sFilename_mosart_domain_out = 'mosart_sag_domain.nc'
    else:
        if sRegion == 'susquehanna':
            sFilename_json_in='/compyfs/liao313/04model/pyhexwatershed/susquehanna/pyhexwatershed20221115001/hexwatershed/hexwatershed.json'
            #sFilename_mpas_in='/qfs/people/liao313/workspace/python/pyhexwatershed_icom/data/susquehanna/input/lnd_cull_mesh.nc'
            sFilename_mosart_parameter_in = '/compyfs/inputdata/rof/mosart/MOSART_Global_half_20210616.nc'
            sFilename_mosart_parameter_out = 'mosart_susquehanna_parameter.nc'
            sFilename_mosart_domain_out = 'mosart_susquehanna_domain.nc'
        else:
            sFilename_json_in='/compyfs/liao313/04model/pyhexwatershed/columbia/pyhexwatershed20221115003/hexwatershed/hexwatershed.json'
            #sFilename_mpas_in='/compyfs/liao313/00raw/mesh/global/lnd_mesh.nc'
            sFilename_mosart_parameter_in = '/compyfs/inputdata/rof/mosart/MOSART_Global_half_20210616.nc'
            sFilename_mosart_parameter_out = 'mosart_columbia_parameter.nc'
            sFilename_mosart_domain_out = 'mosart_columbia_domain.nc'
    
    convert_hexwatershed_json_to_mosart_netcdf(sFilename_json_in, \
        #sFilename_mpas_in, \
            sFilename_mosart_parameter_in,
            sFilename_mosart_parameter_out,\
            sFilename_mosart_domain_out)