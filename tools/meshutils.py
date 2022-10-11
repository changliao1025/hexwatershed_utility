

import json
import numpy as np

# -----------------------------------------------------
# -----------------------------------------------------
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

    with open(MeshJSONfile) as json_file:
        Mesh = json.load(json_file)

    # init lists
    ncell = len(Mesh)
    ID = np.arange(ncell) + 1 # start ID at 1 not 0
    dnID=list()
    cellID = list()
    cellID_downslope=list()

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

# -----------------------------------------------------
# -----------------------------------------------------
def findDownstreamCells(ID,dnID,ipoints,IDtype='mosart'):
    
    """
    find indices and IDs of all cells downstream of each point in a list of points

    Parameters
    ----------
    ID : int (list)
        list of cell IDs (assumes they begin at 1, not 0)
    dnID : int (list)
        list of downstream cell IDs for each element of ID
    ipoints : int (list) | logical (list)
        list of integers representing the indices on ID/dnID for the points of interest
        list of logical true/false. the size (shape) of ipoints must equal the size (shape) of ID and dnID, and is true for the indices of the requested points. 
    IDtype : str
        specifies which cell ID->dnID type. use 'hexwatershed' if ID/dnID are the lCellID and lCellID_downslope fields from the hexwatershed json file. use 'mosart' if ID/dnID are the ID and dnID fields in the mosart parameter file. default: 'mosart'

    Returns
    -------
    i_downstream : int (list)
        list of integers (sublists) representing the indices of all cells downstream of each point in ipoints

    ID_downstream : int (list)
        list of integers (sublists) representing the cell IDs of all cells downstream of each point in ipoints
    
    Author: Matt Cooper (matt.cooper@pnnl.gov)
    """

    ncells = len(ipoints)

    # locate the outlet cell from the global dn_ID
    if IDtype == 'mosart':
        ioutlet = int(np.where(dnID==-9999)[0])
    elif IDtype == 'hexwatershed':
        ioutlet = int(np.where(dnID==-1)[0])

    # init the outputs
    ID_downstream = []
    i_downstream = []

    # loop over all points and find all downstream cells
    for n in range(ncells):
        # idn = ipoints[n] is the cell index that contains this point. ioutlet is the cell index of the outlet. trace down stream from each point until we hit the outlet, collecting the downstream IDs along the way
        idn = ipoints[n]
        dnID_n = []
        dnidx_n = []
        while idn != ioutlet:
            idn = int(np.where(ID==dnID[idn])[0])
            dnID_n.append(int(dnID[idn]))
            dnidx_n.append(idn)

            # note: if ID goes from 0->ncells-1, then use this in place of the assignment above:
            # idn = int(np.where(ID==dnID[idn])[0]) - 1 

        # append these downstream cells to the list for this dam 
        ID_downstream.append(dnID_n)
        i_downstream.append(dnidx_n)

    return i_downstream,ID_downstream
