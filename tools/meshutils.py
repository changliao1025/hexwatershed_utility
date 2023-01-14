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

# -----------------------------------------------------
# -----------------------------------------------------
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



def make_dam_dependency(Dams,Mesh,searchradius,globalID,Mask,FlowLine):

    """
    find dependent cell IDs for each dam (downstream mesh cells that depend on each dam)

    Parameters
    ----------
    Dams : GeoDataFrame
        Points representing computational dams.
    Mesh : GeoDataFrame
        Polygons representing mesh created by pyhexwatershed.
    searchradius : int (scalar)
        horizontal search radius for KDTree algorithm (distance from flowline for dependent cell search)
    globalID : int (list)
        list of cell IDs (assumes they begin at 1, not 0)
    Mask: GeoDataFrame
        LineString (should work if Polygon) representing basin boundary, used to omit mesh cells from dependent cell search
    FlowLine : GeoDataFrame
        Vector flowline. Not required for the algorithm, need to make optional.

    Returns
    -------
    ID_DependentCells : int (list)
        list of integers (sublists) representing the cell IDs of all downstream dependent cells for each dam in Dams
    
    TODO: 
    IDtype : str
        specifies which cell ID->dnID type. use 'hexwatershed' if ID/dnID are the lCellID and lCellID_downslope fields from the hexwatershed json file. use 'mosart' if ID/dnID are the ID and dnID fields in the mosart parameter file. default: 'mosart'

    ipoints : int (list) | logical (list)
        list of integers representing the indices on ID/dnID for the points of interest
        list of logical true/false. the size (shape) of ipoints must equal the size (shape) of ID and dnID, and is true for the indices of the requested points. 
    
    (ipoints would replace Dams)
    
    Author: Matt Cooper (matt.cooper@pnnl.gov), 2023
    """

    import numpy as np
    import inpoly.inpoly2 as inpoly
    from scipy.spatial import KDTree
    try:
        from itertools import izip_longest as zip_longest 
    except ImportError:
        from itertools import zip_longest

    # Part 1: build the kdtree
    # ------------------------

    # add x/y centroid values to Mesh, for the kdtree
    Mesh['X'] = Mesh.centroid.x
    Mesh['Y'] = Mesh.centroid.y
    Dams['X'] = Dams.geometry.x
    Dams['Y'] = Dams.geometry.y

    # build a kdtree for the mesh
    meshTree = KDTree(Mesh[['X','Y']])

    # Part 2: find mesh cells in the basin boundary
    # ---------------------------------------------
    xymesh = np.array(Mesh[['X','Y']]) # mesh cell centroids
    xymask = gdf_coordinate_list(Mask,flatten=False)[0] # boundary
    inmask = np.where(inpoly(xymesh,xymask)[0])[0] # indices
    # inmask = list(inpoly.inpoly2(xymesh,xymask)[0]) # boolean

    # convert to global ID:
    inmask = np.array(globalID.iloc[inmask])

    # Part 3: find mesh cells that contain dams
    # -----------------------------------------

    # find the mesh cells that contain a dam by finding the Mesh (or MeshLine) cells nearest each dam. these are the starting points for the downstream walk to find the dependent cells for each dam
    useflowline = False
    if useflowline is True:
        MeshLine = find_cells_on_vector_flowline(FlowLine,Mesh,meshTree)
        iflowlinedams = KDTree(MeshLine[['X','Y']]).query(Dams[['X','Y']])[1]
        # transform iflowlinedams to the global mesh indices:
        imeshdams = np.array(MeshLine['ID'].iloc[iflowlinedams])
    else:
        imeshdams = meshTree.query(Dams[['X','Y']])[1]

    # add the mesh cell info to the Dams gdf
    Dams['iMesh'] = imeshdams

    # Part 4: find cells downstream of each dam
    # -----------------------------------------

    Dams['i_DownstreamCells'],Dams['ID_DownstreamCells'] = find_downstream_cells(
        np.array(Mesh['ID']),np.array(Mesh['dnID']),imeshdams
        )

    # Part 5: final prep for dependent cell search
    # --------------------------------------------

    # get the mesh cell and dam elevations
    zdams = np.array(Mesh['Elevation'].iloc[imeshdams])
    zmesh = np.array(Mesh[['Elevation']]).flatten()

    # Initialize the DependentCells. Two options: 
    # 1) add an empty ID_DependentCells field to the Dams gdf, or 
    # 2) store in a list (add to gdf outside this algo if necessary)

    # Add to the gdf:
    # Dams['ID_DependentCells'] = np.nan
    # Dams['ID_DependentCells'] = Dams['ID_DependentCells'].astype(object)
    
    # Store in a list:
    ID_DependentCells = []

    # -----------------
    # run the algorithm
    # -----------------

    # starting at each dam, query the tree to find all mesh cells within rxy distance of the flowline below the dam
    for idam,thisdam in Dams.iterrows():

        # use the elevation of the mesh cell that contains the dam as an elevation threshold
        elev = zdams[idam]

        # gather downstream cells for this dam. these are the query points for the kdtree, which finds all cells within rxy distance of these points 
        iquery = Dams['i_DownstreamCells'].iloc[idam]
        MeshQuery = Mesh.iloc[iquery]

        # find all cells within rxy distance of each cell in MeshQuery
        inearby = []
        for icell,thiscell in MeshQuery.iterrows():
            xyquery = (thiscell['X'], thiscell['Y'])
            idx = meshTree.query_ball_point(xyquery,searchradius)
            inearby.extend(idx)

        # keep unique indices that are below elev and in bounds, and convert to globalID
        inearby = np.unique(inearby)
        inearby = globalID.iloc[inearby[np.where(zmesh[inearby] < elev)[0]]]
        inearby = np.intersect1d(inearby,inmask)
        
        # save the cells for this dam
        ID_DependentCells.append(inearby)

        # this puts the Dependent Cells into the data frame (note: use iat to set one element of the gdf not iloc, otherwise we get the slice error)
        # Dams['ID_DependentCells'].iat[idam] = inearby

    # return a uniform-sized padded array
    return np.array(list(zip_longest(*ID_DependentCells, fillvalue=np.nan))).T

    # If adding to the gdf:
    # ID_DependentCells = pd.DataFrame(list(Dams['ID_DependentCells'])).fillna(np.nan).values
    # return ID_DependentCells


def find_cells_on_vector_flowline(FlowLine,Mesh,meshTree):

    """
    find mesh cells that contain a vector flowline vertex

    Parameters
    ----------
    FlowLine : GeoDataFrame
        Vector flowline. Not required for the algorithm, need to make optional.
    Mesh : GeoDataFrame
        Polygons representing mesh created by pyhexwatershed.
    MeshTree : scipy.spatial.kdtree object

    Returns
    -------
    Mesh : GeoDataFrame
        input Mesh GeoDataFrame reduced to cells that contain a flowline vertex (a mesh-based representation of the flowline)
    
    TODO: 
    option to build the kdtree in the function rather than passing it in
    
    Author: Matt Cooper (matt.cooper@pnnl.gov), 2022
    """

    # find mesh cells that contain flowline vertices by finding the nearest mesh cell to each vertex of each line. this won't be necessary if an attribute in Mesh or Line indicates the mapping between them. note: for this query, we need iterrows() b/c each row of Line is a LineString with multiple vertices. 
    imeshline = []
    for iline,thisline in FlowLine.iterrows():
        dist,idx = meshTree.query(thisline.geometry.coords)
        imeshline.append(idx)

    # add the indices to the Line gdf, the flatten the indices and subset the Mesh cells that contain a line segment
    FlowLine['iMesh'] = imeshline
    imeshline = [idx for sublist in imeshline for idx in sublist]

    return Mesh.iloc[imeshline]

    # # for reference, in one line (still needs to be flattened later):
    # # imeshline = [meshTree.query(thisline.geometry.coords)[1] for iline,thisline in Line.iterrows()]

    # Note: this function probably will not ever be necessary. it is only here because the iSegment field in the hexwatershed.json file i used to prototype this code doesn't match the conceptual flowline file (possibly due to my mistake, but I wasn't able to piece it together). we might also keep this if we want the option to build a gdf (called MeshLine in make_dam_dependency) that represents the Mesh cells that contain a flowline, and build a kdtree from that rather than the entire mesh (or an option to identify the flowline from the mesh independently of the iSegment field)


def gdf_coordinate_list(gdf,flatten=False):
    
    '''
    extract coordinate pairs from each feature in a gdf and concatenate them into one list, or list of lists. flatten it into one list if requested.
    Inputs
    gdf : GeoDataFrame
        a geodataframe with features for which the coordinates will be converted to a list
    flatten : boolean
        if true, will return one flattened list of coordinate pairs

    Author: Matt Cooper (matt.cooper@pnnl.gov), 2022
    '''
    from numpy import array
    from shapely.geometry import Polygon

    coordinatelist = []
    for idx,feature in gdf.iterrows():
        
        if all(gdf.geom_type=='LineString'):
            coords = array(feature.geometry.coords)
        elif all(gdf.geom_type=='Point'):
            # not sure if this will work
            coords = array(feature.geometry.coords)
        elif all(gdf.geom_type=='Polygon'):
            coords = array(Polygon(feature['geometry']).exterior.coords)
        
        coordinatelist.append(coords)

        # need a method to figure out if the feature needs to be converted to a polygon first:
        # coords = np.array(Polygon(feature['geometry']).exterior.coords)

    if flatten is True:
        coordinatelist = [coords for features in coordinatelist for coords in features]

    return coordinatelist