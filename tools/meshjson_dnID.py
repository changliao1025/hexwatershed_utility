
def meshjson_dnID(MeshJSONfile):

    """
    find the ID->dnID connectivity from the hexwatershed.json file

    Parameters
    ----------
    MeshJSONfile : str
        Path to the hexwhatershed json file.

    Returns
    -------
    ID : int
        cell ID from 1->ncells
    dnID : int
        downstream cell ID for each ID
    
    Author: Matt Cooper (matt.cooper@pnnl.gov), Donghui Xu and Chang Liao, PNNL
    """
    
    import json
    import numpy as np


    with open(MeshJSONfile) as json_file:
        Mesh = json.load(json_file)

    # init lists
    ncell = len(Mesh)
    ID = np.arange(ncell) + 1 # start ID at 1 not 0
    dnID = list()
    cellID = list()
    cellID_downslope = list()

    for n in range(ncell):
        pcell = Mesh[n]
        cellID.append(int(pcell['lCellID']))
        cellID_downslope.append(int(pcell['lCellID_downslope']))

    #convert to numpy array
    ncell = len(cellID)
    cellID = np.array(cellID)
    cellID_downslope=np.array(cellID_downslope)
    
    for n in range(ncell):
        if cellID_downslope[n] == -1:
            dnID.append(-9999)
        else:
            index = int(np.where(cellID == cellID_downslope[n])[0])
            dnID.append(index + 1)
        
    dnID = np.array(dnID)
    return ID,dnID

# started to make one that works with the Mesh geodataframe, but it isnt' needed for now
# def meshgdf_dnID(Mesh):
    # for icell,thiscell in Mesh.iterrows():
