import os
from osgeo import ogr

def get_outlet_location(sFilename_simplified_flowline_in):
    """
    Get the outlet location (downstream endpoint) from a river network file.

    Args:
        sFilename_simplified_flowline_in (str): Path to the input vector file (GeoJSON, shapefile, etc.)
        iFlag_furthest_downstream (int): If 1, find the furthest downstream point across all features.
                                         If 0, return the last point of the first feature.

    Returns:
        tuple: (longitude, latitude) of the outlet location

    Raises:
        ValueError: If file doesn't exist, cannot be opened, or has no valid geometry
    """
    # Check if the file exists
    if not os.path.exists(sFilename_simplified_flowline_in):
        raise ValueError(f"File {sFilename_simplified_flowline_in} does not exist")

    # Open the dataset
    pDataset = ogr.Open(sFilename_simplified_flowline_in)
    if pDataset is None:
        raise ValueError(f"Could not open {sFilename_simplified_flowline_in}")

    try:
        pLayer = pDataset.GetLayer(0)
        if pLayer is None:
            raise ValueError(f"Could not get layer from {sFilename_simplified_flowline_in}")

        nFeature_count = pLayer.GetFeatureCount()
        if nFeature_count == 0:
            raise ValueError(f"No features found in {sFilename_simplified_flowline_in}")


        # Simple case: return last point of first feature
        pFeature = pLayer.GetFeature(0)
        if pFeature is None:
            raise ValueError(f"Could not get first feature from {sFilename_simplified_flowline_in}")
        pGeometry = pFeature.GetGeometryRef()
        if pGeometry is None or pGeometry.GetPointCount() == 0:
            raise ValueError(f"First feature has no valid geometry")
        # Get the last point of the geometry (outlet)
        pPoint = pGeometry.GetPoint(pGeometry.GetPointCount() - 1)
        dLongitude = pPoint[0]
        dLatitude = pPoint[1]

        return dLongitude, dLatitude

    finally:
        # Clean up
        pDataset = None