import logging
import os
from typing import Dict, Optional

from osgeo import ogr, osr
from pyearth.gis.gdal.gdal_vector_format_support import get_vector_driver_from_filename
from pyearth.toolbox.mesh.vertex import pyvertex
from pyflowline.formats.read_flowline import read_flowline_vector
from rtree.index import Index as RTreeindex

logger = logging.getLogger(__name__)


def snap_dams_to_river_network(
    sFilename_dam: str,
    sFilename_river: str,
    sFilename_out: str,
    sVariable_drainage: Optional[str] = "drainage",
    sVariable_id: Optional[str] = None,
    dDistance_max_snap: float = 3000.0,
    dDistance_on_river: float = 10.0,
    dBuffer_search_deg: float = 0.05,
) -> Dict[str, str]:
    """
    Snap dam locations to the nearest vertex of a river flowline network.

    Categorizes dams into:
    1. 'on': Dam is already on the river network (distance <= dDistance_on_river).
    2. 'snapped': Dam is successfully snapped to nearest river vertex (distance <= dDistance_max_snap).
    3. 'no_river': No river candidate within initial bounding box search radius.
    4. 'too_far': Nearest river flowline exceeds dDistance_max_snap.
    5. 'merged': Combined 'on' and 'snapped' dams on the network.

    Parameters
    ----------
    sFilename_dam : str
        Path to input dam point vector file (Shapefile, GeoJSON, etc.).
    sFilename_river : str
        Path to input river network flowline file.
    sFilename_out : str
        Base path for output files.
    sVariable_drainage : Optional[str], default='drainage'
        Attribute field name containing drainage area. If None or not found, set to 0.0.
    sVariable_id : Optional[str], default=None
        Optional field name for dam ID. If None, auto-enumerated.
    dDistance_max_snap : float, default=3000.0
        Maximum allowable snapping distance in meters.
    dDistance_on_river : float, default=10.0
        Threshold distance in meters to consider dam already situated on river.
    dBuffer_search_deg : float, default=0.05
        Spatial bounding-box search radius in degrees (~5 km).

    Returns
    -------
    Dict[str, str]
        Dictionary mapping category keys ('on', 'snapped', 'merged', 'no_river', 'too_far')
        to their respective output file paths.
    """
    if not os.path.exists(sFilename_dam):
        raise FileNotFoundError(f"Dam location file does not exist: {sFilename_dam}")

    if not os.path.exists(sFilename_river):
        raise FileNotFoundError(f"River network file does not exist: {sFilename_river}")

    sDir_out = os.path.dirname(sFilename_out)
    if sDir_out and not os.path.exists(sDir_out):
        os.makedirs(sDir_out, exist_ok=True)

    pDataSource_dam = ogr.Open(sFilename_dam, 0)
    if pDataSource_dam is None:
        raise RuntimeError(f"Could not open dam dataset: {sFilename_dam}")
    pLayer_dam = pDataSource_dam.GetLayer()
    nDam = pLayer_dam.GetFeatureCount()

    srs_target = osr.SpatialReference()
    srs_target.ImportFromEPSG(4326)

    # Output filenames for scenarios
    sFilename_base, sExt = os.path.splitext(sFilename_out)
    sFilename_on = f"{sFilename_base}_on{sExt}"
    sFilename_snapped = f"{sFilename_base}_snapped{sExt}"
    sFilename_merged = f"{sFilename_base}_merged{sExt}"
    sFilename_no_river = f"{sFilename_base}_no_river{sExt}"
    sFilename_too_far = f"{sFilename_base}_too_far{sExt}"

    pDriver = get_vector_driver_from_filename(sFilename_out)
    if pDriver is None:
        raise ValueError(f"Could not determine driver for: {sFilename_out}")

    for sFile in [
        sFilename_on,
        sFilename_snapped,
        sFilename_merged,
        sFilename_no_river,
        sFilename_too_far,
    ]:
        if os.path.exists(sFile):
            pDriver.DeleteDataSource(sFile)

    pDataSource_on = pDriver.CreateDataSource(sFilename_on)
    pLayer_on = pDataSource_on.CreateLayer("dam", srs_target, ogr.wkbPoint)

    pDataSource_snapped = pDriver.CreateDataSource(sFilename_snapped)
    pLayer_snapped = pDataSource_snapped.CreateLayer("dam", srs_target, ogr.wkbPoint)

    pDataSource_merged = pDriver.CreateDataSource(sFilename_merged)
    pLayer_merged = pDataSource_merged.CreateLayer("dam", srs_target, ogr.wkbPoint)

    pDataSource_no_river = pDriver.CreateDataSource(sFilename_no_river)
    pLayer_no_river = pDataSource_no_river.CreateLayer("dam", srs_target, ogr.wkbPoint)

    pDataSource_too_far = pDriver.CreateDataSource(sFilename_too_far)
    pLayer_too_far = pDataSource_too_far.CreateLayer("dam", srs_target, ogr.wkbPoint)

    for pLayer in [pLayer_snapped, pLayer_no_river, pLayer_too_far, pLayer_on, pLayer_merged]:
        pLayer.CreateField(ogr.FieldDefn("id", ogr.OFTInteger))
        pLayer.CreateField(ogr.FieldDefn("dLongitude_degree", ogr.OFTReal))
        pLayer.CreateField(ogr.FieldDefn("dLatitude_degree", ogr.OFTReal))
        pLayer.CreateField(ogr.FieldDefn("drainage", ogr.OFTReal))

    pLayerDefn_snapped = pLayer_snapped.GetLayerDefn()
    pLayerDefn_no_river = pLayer_no_river.GetLayerDefn()
    pLayerDefn_too_far = pLayer_too_far.GetLayerDefn()
    pLayerDefn_on = pLayer_on.GetLayerDefn()
    pLayerDefn_merged = pLayer_merged.GetLayerDefn()

    aFlowline, _ = read_flowline_vector(sFilename_river)
    nFlowline = len(aFlowline)

    index_flowline = RTreeindex()
    for i in range(nFlowline):
        pBound = aFlowline[i].pBound
        index_flowline.insert(i, pBound)

    # Check dam attribute field names
    dam_defn = pLayer_dam.GetLayerDefn()
    dam_field_names = [dam_defn.GetFieldDefn(j).GetName() for j in range(dam_defn.GetFieldCount())]

    has_drainage = sVariable_drainage in dam_field_names if sVariable_drainage else False
    has_id = sVariable_id in dam_field_names if sVariable_id else False

    iDam_on = 1
    iDam_snapped = 1
    iDam_no_river = 1
    iDam_too_far = 1
    nCount_on = 0
    nCount_snapped = 0
    nCount_no_river = 0
    nCount_too_far = 0

    pLayer_dam.ResetReading()
    for i, pFeature_dam in enumerate(pLayer_dam):
        point = pFeature_dam.GetGeometryRef()
        if point is None:
            continue

        dDrainage = float(pFeature_dam.GetField(sVariable_drainage)) if has_drainage else 0.0
        nDam_id = int(pFeature_dam.GetField(sVariable_id)) if has_id else i + 1

        x = point.GetX()
        y = point.GetY()
        pVertex = pyvertex({"dLongitude_degree": x, "dLatitude_degree": y})

        pBound = (
            x - dBuffer_search_deg,
            y - dBuffer_search_deg,
            x + dBuffer_search_deg,
            y + dBuffer_search_deg,
        )
        aIntersect = list(index_flowline.intersection(pBound))

        if len(aIntersect) == 0:
            # Scenario: Dam is not close to any candidate flowline
            pPoint = ogr.Geometry(ogr.wkbPoint)
            pPoint.AddPoint(x, y)
            pFeature_out = ogr.Feature(pLayerDefn_no_river)
            pFeature_out.SetGeometry(pPoint)
            pFeature_out.SetField("id", nDam_id)
            pFeature_out.SetField("dLongitude_degree", x)
            pFeature_out.SetField("dLatitude_degree", y)
            pFeature_out.SetField("drainage", dDrainage)
            pLayer_no_river.CreateFeature(pFeature_out)
            iDam_no_river += 1
            nCount_no_river += 1
        else:
            dDistance_min = 1.0e10
            pVertex_closest = None

            for k in aIntersect:
                pFlowline = aFlowline[k]
                distance, pVertex_out = pFlowline.calculate_distance_to_vertex(pVertex)
                if distance < dDistance_min:
                    dDistance_min = distance
                    pVertex_closest = pVertex_out

            if dDistance_min <= dDistance_max_snap and pVertex_closest is not None:
                snap_lon = pVertex_closest.dLongitude_degree
                snap_lat = pVertex_closest.dLatitude_degree

                if dDistance_min <= dDistance_on_river:
                    # Already on river network
                    target_layer = pLayer_on
                    target_defn = pLayerDefn_on
                    nCount_on += 1
                    iDam_on += 1
                else:
                    # Snapped within tolerance
                    target_layer = pLayer_snapped
                    target_defn = pLayerDefn_snapped
                    nCount_snapped += 1
                    iDam_snapped += 1

                pPoint = ogr.Geometry(ogr.wkbPoint)
                pPoint.AddPoint(snap_lon, snap_lat)
                pFeature_out = ogr.Feature(target_defn)
                pFeature_out.SetGeometry(pPoint)
                pFeature_out.SetField("id", nDam_id)
                pFeature_out.SetField("dLongitude_degree", snap_lon)
                pFeature_out.SetField("dLatitude_degree", snap_lat)
                pFeature_out.SetField("drainage", dDrainage)
                target_layer.CreateFeature(pFeature_out)

                # Add to merged on-river dataset
                pFeature_merged = ogr.Feature(pLayerDefn_merged)
                pFeature_merged.SetGeometry(pPoint)
                pFeature_merged.SetField("id", nDam_id)
                pFeature_merged.SetField("dLongitude_degree", snap_lon)
                pFeature_merged.SetField("dLatitude_degree", snap_lat)
                pFeature_merged.SetField("drainage", dDrainage)
                pLayer_merged.CreateFeature(pFeature_merged)

            else:
                # Exceeds max snap distance
                pPoint = ogr.Geometry(ogr.wkbPoint)
                pPoint.AddPoint(x, y)
                pFeature_out = ogr.Feature(pLayerDefn_too_far)
                pFeature_out.SetGeometry(pPoint)
                pFeature_out.SetField("id", nDam_id)
                pFeature_out.SetField("dLongitude_degree", x)
                pFeature_out.SetField("dLatitude_degree", y)
                pFeature_out.SetField("drainage", dDrainage)
                pLayer_too_far.CreateFeature(pFeature_out)
                iDam_too_far += 1
                nCount_too_far += 1

    pDataSource_dam = None
    pDataSource_snapped = None
    pDataSource_no_river = None
    pDataSource_on = None
    pDataSource_too_far = None
    pDataSource_merged = None

    logger.info("Total dams processed: %d", nDam)
    logger.info(
        "Already on river (<= %.1fm): %d -> %s", dDistance_on_river, nCount_on, sFilename_on
    )
    logger.info(
        "Successfully snapped (<= %.1fm): %d -> %s",
        dDistance_max_snap,
        nCount_snapped,
        sFilename_snapped,
    )
    logger.info(
        "Too far (> %.1fm): %d -> %s", dDistance_max_snap, nCount_too_far, sFilename_too_far
    )
    logger.info("No river candidates: %d -> %s", nCount_no_river, sFilename_no_river)

    return {
        "on": sFilename_on,
        "snapped": sFilename_snapped,
        "merged": sFilename_merged,
        "no_river": sFilename_no_river,
        "too_far": sFilename_too_far,
    }


if __name__ == "__main__":
    import sys

    if len(sys.argv) > 3:
        snap_dams_to_river_network(sys.argv[1], sys.argv[2], sys.argv[3])
