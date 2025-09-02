
def find_downstream_cells(ID,dnID,ipoints,IDtype='mosart'):
    
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
    
    # imports
    import numpy as np

    # function definitions
    intvector = np.vectorize(np.int_)

    # initialize
    ncells = len(ipoints)
    ID_downstream = []
    i_downstream = []

    # locate the outlet cell from the global dn_ID
    if IDtype == 'mosart':
        ioutlet = intvector(np.where(dnID==-9999)[0])
    elif IDtype == 'hexwatershed':
        ioutlet = intvector(np.where(dnID==-1)[0])

    # loop over all points and find all downstream cells
    for n in range(ncells):
        # idn = ipoints[n] is the cell index that contains this point. ioutlet is the cell index of the outlet. trace down stream from each point until we hit the outlet, collecting the downstream IDs along the way
        idn = ipoints[n]
        dnID_n = []
        dnidx_n = []
        while idn not in ioutlet:  # while idn != ioutlet: (if ioutlet is scalar)
            idn = int(np.where(ID==dnID[idn])[0])
            dnID_n.append(int(dnID[idn]))
            dnidx_n.append(idn)

            # note: if ID goes from 0->ncells-1, then use this in place of the assignment above:
            # idn = int(np.where(ID==dnID[idn])[0]) - 1 

        # append these downstream cells to the list for this dam 
        ID_downstream.append(dnID_n)
        i_downstream.append(dnidx_n)

    return i_downstream,ID_downstream
