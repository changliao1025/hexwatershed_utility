import os
from typing import Optional, Tuple

from osgeo import ogr


def get_outlet_location(
    sFilename_simplified_flowline_in: str,
    iFlag_furthest_downstream: int = 0,
    sVariable_downstream: Optional[str] = "downstream_id",
) -> Tuple[float, float]:
    """
    Get the outlet location (longitude, latitude) from a river network vector file.

    Parameters
    ----------
    sFilename_simplified_flowline_in : str
        Path to input river network vector file (GeoJSON, Shapefile, etc.).
    iFlag_furthest_downstream : int, default=0
        If 1, identify the furthest downstream terminal reach (outlet) across all features.
        If 0, return the last vertex of the first feature in the file.
    sVariable_downstream : Optional[str], default='downstream_id'
        Field name representing downstream reach ID used when iFlag_furthest_downstream=1.

    Returns
    -------
    Tuple[float, float]
        (longitude, latitude) of the river network outlet.

    Raises
    ------
    FileNotFoundError
        If input file does not exist.
    ValueError
        If dataset cannot be opened or contains no valid geometries.
    """
    if not os.path.exists(sFilename_simplified_flowline_in):
        raise FileNotFoundError(f"File {sFilename_simplified_flowline_in} does not exist")

    pDataset = ogr.Open(sFilename_simplified_flowline_in, 0)
    if pDataset is None:
        raise ValueError(f"Could not open {sFilename_simplified_flowline_in}")

    try:
        pLayer = pDataset.GetLayer(0)
        if pLayer is None:
            raise ValueError(f"Could not get layer from {sFilename_simplified_flowline_in}")

        nFeature_count = pLayer.GetFeatureCount()
        if nFeature_count == 0:
            raise ValueError(f"No features found in {sFilename_simplified_flowline_in}")

        pTargetFeature = None

        if iFlag_furthest_downstream == 1:
            # Search for terminal reach (e.g. downstream_id <= 0 or not found)
            pLayerDefn = pLayer.GetLayerDefn()
            field_names = [
                pLayerDefn.GetFieldDefn(i).GetName() for i in range(pLayerDefn.GetFieldCount())
            ]

            dn_col = sVariable_downstream if sVariable_downstream in field_names else None
            pLayer.ResetReading()
            for pFeature in pLayer:
                if dn_col is not None:
                    val = pFeature.GetField(dn_col)
                    if val is None or val in (0, -1, "-1", "0", ""):
                        pTargetFeature = pFeature
                        break

            if pTargetFeature is None:
                # Fallback to feature with max FID or last feature
                pTargetFeature = pLayer.GetFeature(nFeature_count - 1)
        else:
            pTargetFeature = pLayer.GetFeature(0)

        if pTargetFeature is None:
            raise ValueError(f"Could not retrieve feature from {sFilename_simplified_flowline_in}")

        pGeometry = pTargetFeature.GetGeometryRef()
        if pGeometry is None or pGeometry.GetPointCount() == 0:
            raise ValueError("Feature has no valid geometry points")

        pPoint = pGeometry.GetPoint(pGeometry.GetPointCount() - 1)
        dLongitude = float(pPoint[0])
        dLatitude = float(pPoint[1])

        return dLongitude, dLatitude

    finally:
        pDataset = None


__all__ = [
    "get_outlet_location",
]
