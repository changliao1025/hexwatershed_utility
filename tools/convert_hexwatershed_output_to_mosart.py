import json
def convert_hexwatershed_json_to_mosart_netcdf(sFilename_json_in, sFilename_mpas_in, sFilename_mosart_in):
    #open json and read data
    aCellID = list()
    aCellID_downslope=list()
    with open(sFilename_json_in) as json_file:
        data = json.load(json_file)      
        ncell = len(data)
        lID =0 
        for i in range(ncell):
            pcell = data[i]
            
            aCellID.append(int(pcell['lCellID']))
            aCellID_downslope.append(int(pcell['lCellID_downslope']))
            aElevation.append(float(pcell['dElevation']))
            aDrainage.append(float(pcell['dElevation']))
        pass
    
    #open mpas netcdf 
    
    

    return
def prepare_mosart_domain_file():
    return 