import os
from osgeo import gdal, ogr
gdal.UseExceptions()
import importlib.util
iFlag_cython = importlib.util.find_spec("cython")
if iFlag_cython is not None:
    from tinyr import RTree
    iFlag_use_rtree = 1
else:
    iFlag_use_rtree = 0

from pyearth.toolbox.geometry.create_gcs_buffer_zone import create_polyline_buffer_zone
from pyearth.gis.geometry.calculate_distance_based_on_longitude_latitude import calculate_distance_based_on_longitude_latitude
from pyearth.toolbox.management.vector.merge_features import merge_features

def geometries_bbox_overlap(bbox1, bbox2, tolerance=1e-10):
    """
    Check if two bounding boxes overlap with optional tolerance

    :param bbox1: (minX, maxX, minY, maxY) for first geometry
    :param bbox2: (minX, maxX, minY, maxY) for second geometry
    :param tolerance: small buffer for near-touching geometries
    :return: True if bounding boxes overlap
    """
    return not (bbox1[1] + tolerance < bbox2[0] or  # bbox1.maxX < bbox2.minX
                bbox2[1] + tolerance < bbox1[0] or  # bbox2.maxX < bbox1.minX
                bbox1[3] + tolerance < bbox2[2] or  # bbox1.maxY < bbox2.minY
                bbox2[3] + tolerance < bbox1[2])    # bbox2.maxY < bbox1.minY

def fix_naturalearth_hydrosheds_incompatibility(aFilename_hydrosheds_flowline, sFilename_vector_naturalearth, sFilename_vector_naturalearth_updated):
    """
    Fix incompatibility between Natural Earth land polygons and HydroSHEDS flowlines.

    This function identifies flowlines that are outside the Natural Earth land polygon
    and creates buffer zones around them to extend the land polygon, ensuring consistency
    between the land mask and river network. Uses RTree for efficient spatial indexing.

    Args:
        aFilename_hydrosheds_flowline: List of HydroSHEDS flowline file paths
        sFilename_vector_naturalearth: Path to Natural Earth land polygon file
        sFilename_vector_naturalearth_updated: Path for the updated Natural Earth file

    Returns:
        sFilename_vector_naturalearth_updated: Path to the updated file
    """

    nFlowline = len(aFilename_hydrosheds_flowline)

    # Open the Natural Earth land polygon dataset
    pDataset_naturalearth = ogr.Open(sFilename_vector_naturalearth)
    if pDataset_naturalearth is None:
        print(f'Error: Could not open Natural Earth file {sFilename_vector_naturalearth}')
        return None

    pLayer_naturalearth = pDataset_naturalearth.GetLayer(0)
    if pLayer_naturalearth is None:
        print('Error: No layer found in Natural Earth file')
        return None

    # Get spatial reference from Natural Earth layer
    pSpatial_reference = pLayer_naturalearth.GetSpatialRef()

    # Build spatial index for Natural Earth land polygons
    print('Building spatial index for Natural Earth land polygons...')
    land_geometries = []
    land_spatial_index = None

    if iFlag_use_rtree:
        print('Using RTree for spatial indexing')
        land_spatial_index = RTree()
    else:
        print('Using dictionary-based spatial indexing (install tinyr for better performance)')
        land_spatial_index = {}

    # Index all land polygons
    nLand_features = pLayer_naturalearth.GetFeatureCount()
    for idx, pFeature_land in enumerate(pLayer_naturalearth):
        pGeometry_land = pFeature_land.GetGeometryRef()
        if pGeometry_land is None:
            continue

        land_geometries.append(pGeometry_land.Clone())
        envelope = pGeometry_land.GetEnvelope()  # (minX, maxX, minY, maxY)

        if iFlag_use_rtree:
            # RTree expects (minX, minY, maxX, maxY)
            bbox = (envelope[0], envelope[2], envelope[1], envelope[3])
            land_spatial_index.insert(idx, bbox)
        else:
            land_spatial_index[idx] = envelope

    print(f'Indexed {len(land_geometries)} land polygon features')

    # Create a list to store buffer geometries for flowlines outside land
    aBuffer_geometries = []

    for iFlowline in range(nFlowline):
        sFilename_flowline = aFilename_hydrosheds_flowline[iFlowline]

        print(f'Processing flowline file {iFlowline + 1}/{nFlowline}: {os.path.basename(sFilename_flowline)}')

        # Open the flowline dataset
        pDataset_flowline = ogr.Open(sFilename_flowline)
        if pDataset_flowline is None:
            print(f'Warning: Could not open flowline file {sFilename_flowline}')
            continue

        pLayer_flowline = pDataset_flowline.GetLayer(0)
        if pLayer_flowline is None:
            print(f'Warning: No layer found in flowline file {sFilename_flowline}')
            continue

        # Process each flowline feature
        nFlowline_features = pLayer_flowline.GetFeatureCount()
        nOutside_count = 0

        for pFeature_flowline in pLayer_flowline:
            pGeometry_flowline = pFeature_flowline.GetGeometryRef()
            if pGeometry_flowline is None:
                continue

            # Get flowline bounding box for spatial query
            flowline_envelope = pGeometry_flowline.GetEnvelope()  # (minX, maxX, minY, maxY)

            # Find candidate land polygons using spatial index
            candidate_indices = []

            if iFlag_use_rtree:
                # RTree expects (minX, minY, maxX, maxY)
                query_bbox = (flowline_envelope[0], flowline_envelope[2], flowline_envelope[1], flowline_envelope[3])
                candidate_indices = list(land_spatial_index.search(query_bbox))
            else:
                # Fallback: check bounding box overlap for all polygons
                for idx, land_envelope in land_spatial_index.items():
                    if geometries_bbox_overlap(flowline_envelope, land_envelope):
                        candidate_indices.append(idx)

            # Check if flowline intersects with any candidate land polygon
            bIntersects_land = False

            # Start with the original flowline
            whole_flowline = pGeometry_flowline.Clone()

            for candidate_idx in candidate_indices:
                if candidate_idx >= len(land_geometries):
                    continue

                pGeometry_land = land_geometries[candidate_idx]
                if pGeometry_land is None:
                    continue


                # Check if flowline intersects with land polygon
                if whole_flowline.Intersects(pGeometry_land):
                    bIntersects_land = True
                    # Get the part of flowline that is outside this land polygon
                    outside_part = whole_flowline.Difference(pGeometry_land)
                    if outside_part is not None and not outside_part.IsEmpty():
                        # Update remaining flowline to only the outside part
                        outside_flowline = outside_part.Clone()
                        # Check if the outside part has significant length
                        #use create_polyline_buffer_zone to get the buffer
                        outside_flowline.FlattenTo2D()
                        sGeometry_type = outside_flowline.GetGeometryName()
                        if sGeometry_type == 'LINESTRING':
                            nOutside_count += 1
                            sWkt= outside_flowline.ExportToWkt()
                            #calcuate the length of the outside_part to determine the buffer distance
                            point_start = outside_flowline.GetPoint(0)
                            point_end = outside_flowline.GetPoint(outside_flowline.GetPointCount() - 1)
                            distance = calculate_distance_based_on_longitude_latitude(point_start[0], point_start[1], point_end[0], point_end[1])
                            sWkt_buffer_polygon = create_polyline_buffer_zone(sWkt, distance/2.0)
                            #create a geometry from the buffer wkt
                            pBuffer_geometry = ogr.CreateGeometryFromWkt(sWkt_buffer_polygon)
                            aBuffer_geometries.append(pBuffer_geometry)
                        else:
                            if sGeometry_type == 'MULTILINESTRING':
                                for i in range(outside_flowline.GetGeometryCount()):
                                    outside_part_i = outside_flowline.GetGeometryRef(i)
                                    if outside_part_i is None or outside_part_i.IsEmpty():
                                        continue
                                    nOutside_count += 1
                                    sWkt= outside_part_i.ExportToWkt()
                                    point_start = outside_part_i.GetPoint(0)
                                    point_end = outside_part_i.GetPoint(outside_part_i.GetPointCount() - 1)
                                    distance = calculate_distance_based_on_longitude_latitude(point_start[0], point_start[1], point_end[0], point_end[1])
                                    #calcuate the length of the outside_part to determine the buffer distance
                                    sWkt_buffer_polygon = create_polyline_buffer_zone(sWkt, distance/2.0)
                                    #create a geometry from the buffer wkt
                                    pBuffer_geometry = ogr.CreateGeometryFromWkt(sWkt_buffer_polygon)
                                    aBuffer_geometries.append(pBuffer_geometry)
                        #use this distance to create a circle geometry as buffer
                        pass

                # If bIntersects_land is True but no outside_flowline_parts,
                # it means the entire flowline is inside land polygons - no buffer needed

        print(f'  Found {nOutside_count} flowline features (or parts) requiring buffers out of {nFlowline_features} total')

        # Clean up flowline dataset
        pDataset_flowline = None

    # Create the updated Natural Earth file

    #first create a tmp file to store the original natural earth features
    sExtension = os.path.splitext(sFilename_vector_naturalearth)[1]
    #get the base name without extension
    sBase_name = os.path.splitext(os.path.basename(sFilename_vector_naturalearth))[0]

    sFilename_vector_naturalearth_tmp = os.path.join(os.path.dirname(sFilename_vector_naturalearth), sBase_name + '_tmp'+ sExtension)


    sFormat = 'GeoJSON'  # Default format

    pDriver = ogr.GetDriverByName(sFormat)
    if pDriver is None:
        print(f'Error: Driver {sFormat} not available')
        return None

    # Remove existing output file if it exists
    if os.path.exists(sFilename_vector_naturalearth_tmp):
        pDriver.DeleteDataSource(sFilename_vector_naturalearth_tmp)

    # Create output dataset
    pDataset_output = pDriver.CreateDataSource(sFilename_vector_naturalearth_tmp)
    if pDataset_output is None:
        print('Error: Could not create output dataset')
        return None

    # Create output layer
    pLayer_output = pDataset_output.CreateLayer('naturalearth_updated', pSpatial_reference, geom_type=ogr.wkbPolygon)
    if pLayer_output is None:
        print('Error: Could not create output layer')
        return None

    # Copy field definitions from input layer
    pLayerDefn_input = pLayer_naturalearth.GetLayerDefn()
    for i in range(pLayerDefn_input.GetFieldCount()):
        pFieldDefn = pLayerDefn_input.GetFieldDefn(i)
        if pLayer_output.CreateField(pFieldDefn) != 0:
            print(f'Warning: Failed to create field {pFieldDefn.GetName()}')

    # Copy original Natural Earth features to output
    pLayer_naturalearth.ResetReading()
    for pFeature in pLayer_naturalearth:
        pFeature_out = ogr.Feature(pLayer_output.GetLayerDefn())
        pFeature_out.SetFrom(pFeature)
        if pLayer_output.CreateFeature(pFeature_out) != 0:
            print('Warning: Failed to create feature in output layer')
        pFeature_out = None


    # Add buffer geometries as new features
    if aBuffer_geometries:
        print(f'Adding {len(aBuffer_geometries)} buffer zones for flowlines outside land')
        for pBuffer_geom in aBuffer_geometries:
            if pBuffer_geom is None or pBuffer_geom.IsEmpty():
                continue

            pFeature_out = ogr.Feature(pLayer_output.GetLayerDefn())
            pFeature_out.SetGeometry(pBuffer_geom)
            if pLayer_output.CreateFeature(pFeature_out) != 0:
                print('Warning: Failed to create buffer feature in output layer')
            pFeature_out = None

    # Clean up
    pDataset_output = None
    pDataset_naturalearth = None

    merge_features(sFilename_vector_naturalearth_tmp, sFilename_vector_naturalearth_updated)

    print(f'Successfully created updated Natural Earth file: {sFilename_vector_naturalearth_updated}')
    print(f'Added {len(aBuffer_geometries)} buffer zones for rivers outside original land polygons')

    return
