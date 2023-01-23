import numpy as np
import inpoly.inpoly2 as inpoly
from scipy.spatial import KDTree
try:
    from itertools import izip_longest as zip_longest 
except ImportError:
    from itertools import zip_longest
from shapely.geometry import Polygon


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

    coordinatelist = []
    for idx,feature in gdf.iterrows():
        
        if all(gdf.geom_type=='LineString'):
            coords = np.array(feature.geometry.coords)
        elif all(gdf.geom_type=='Point'):
            # not sure if this will work
            coords = np.array(feature.geometry.coords)
        elif all(gdf.geom_type=='Polygon'):
            coords = np.array(Polygon(feature['geometry']).exterior.coords)
        
        coordinatelist.append(coords)

        # need a method to figure out if the feature needs to be converted to a polygon first:
        # coords = np.array(Polygon(feature['geometry']).exterior.coords)

    if flatten is True:
        coordinatelist = [coords for features in coordinatelist for coords in features]

    return coordinatelist