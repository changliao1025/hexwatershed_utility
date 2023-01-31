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