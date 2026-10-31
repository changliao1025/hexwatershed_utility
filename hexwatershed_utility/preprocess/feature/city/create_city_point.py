import csv
import logging
import os
from typing import Optional

from osgeo import ogr, osr
from pyearth.gis.gdal.gdal_vector_format_support import get_vector_driver_from_filename

logger = logging.getLogger(__name__)


def _find_column_index(headers, candidate_names, default_idx=None):
    """Find index of a column matching candidate names (case-insensitive)."""
    headers_lower = [h.strip().lower() for h in headers]
    for name in candidate_names:
        if name.lower() in headers_lower:
            return headers_lower.index(name.lower())
    return default_idx


def create_city_point(
    sFilename_in: str,
    sFilename_out: str,
    dPopulation_min: float = 1.0e6,
    sCol_city: Optional[str] = None,
    sCol_lat: Optional[str] = None,
    sCol_lon: Optional[str] = None,
    sCol_pop: Optional[str] = None,
) -> str:
    """
    Extract major cities exceeding a population threshold into a vector point dataset.

    Supports both CSV files (e.g. SimpleMaps World Cities, GeoNames) and vector files
    (e.g. Natural Earth populated places Shapefile/GeoJSON).

    Parameters
    ----------
    sFilename_in : str
        Path to input city dataset (CSV or vector format like Shapefile, GeoJSON, GeoPackage).
    sFilename_out : str
        Path to output vector file (GeoJSON, Shapefile, etc.).
    dPopulation_min : float, default=1.0e6
        Minimum population cutoff.
    sCol_city : Optional[str], default=None
        Column/field name for city name. If None, auto-detected.
    sCol_lat : Optional[str], default=None
        Column/field name for latitude (CSV inputs). If None, auto-detected.
    sCol_lon : Optional[str], default=None
        Column/field name for longitude (CSV inputs). If None, auto-detected.
    sCol_pop : Optional[str], default=None
        Column/field name for population. If None, auto-detected.

    Returns
    -------
    str
        Path to generated output vector file.
    """
    if not os.path.exists(sFilename_in):
        raise FileNotFoundError(f"Input city dataset not found: {sFilename_in}")

    sDir_out = os.path.dirname(sFilename_out)
    if sDir_out and not os.path.exists(sDir_out):
        os.makedirs(sDir_out, exist_ok=True)

    pDriver_out = get_vector_driver_from_filename(sFilename_out)
    if pDriver_out is None:
        raise ValueError(f"Could not get vector driver for output: {sFilename_out}")

    if os.path.exists(sFilename_out):
        pDriver_out.DeleteDataSource(sFilename_out)

    srs_target = osr.SpatialReference()
    srs_target.ImportFromEPSG(4326)

    pDataSource_out = pDriver_out.CreateDataSource(sFilename_out)
    if pDataSource_out is None:
        raise RuntimeError(f"Could not create output dataset: {sFilename_out}")

    pLayer_out = pDataSource_out.CreateLayer("city", srs_target, ogr.wkbPoint)
    pLayer_out.CreateField(ogr.FieldDefn("city", ogr.OFTString))
    pLayer_out.CreateField(ogr.FieldDefn("latitude", ogr.OFTReal))
    pLayer_out.CreateField(ogr.FieldDefn("longitude", ogr.OFTReal))
    pLayer_out.CreateField(ogr.FieldDefn("population", ogr.OFTReal))
    pLayerDefn_out = pLayer_out.GetLayerDefn()

    count = 0

    if sFilename_in.lower().endswith(".csv"):
        with open(sFilename_in, encoding="utf-8", errors="ignore") as f:
            reader = csv.reader(f)
            header = next(reader, None)
            if header is None:
                pDataSource_out = None
                return sFilename_out

            city_idx = (
                _find_column_index(header, [sCol_city])
                if sCol_city
                else _find_column_index(
                    header,
                    ["city", "city_ascii", "name", "city_name", "place_name", "asciiname"],
                    default_idx=0,
                )
            )
            lat_idx = (
                _find_column_index(header, [sCol_lat])
                if sCol_lat
                else _find_column_index(header, ["lat", "latitude", "y", "lat_deg"], default_idx=2)
            )
            lon_idx = (
                _find_column_index(header, [sCol_lon])
                if sCol_lon
                else _find_column_index(
                    header, ["lng", "lon", "long", "longitude", "x", "lon_deg"], default_idx=3
                )
            )
            pop_idx = (
                _find_column_index(header, [sCol_pop])
                if sCol_pop
                else _find_column_index(
                    header, ["population", "pop", "pop_max", "pop_min"], default_idx=9
                )
            )

            max_needed = max(
                idx for idx in (city_idx, lat_idx, lon_idx, pop_idx) if idx is not None
            )

            for row in reader:
                if len(row) <= max_needed:
                    continue

                sName = row[city_idx] if city_idx is not None else ""
                try:
                    dLat = float(row[lat_idx])
                    dLon = float(row[lon_idx])
                    pop_str = row[pop_idx].strip() if pop_idx is not None else "0"
                    dPop = float(pop_str) if pop_str else 0.0
                except (ValueError, IndexError):
                    continue

                if dPop < dPopulation_min:
                    continue

                pPoint = ogr.Geometry(ogr.wkbPoint)
                pPoint.AddPoint(dLon, dLat)

                pFeature = ogr.Feature(pLayerDefn_out)
                pFeature.SetGeometry(pPoint)
                pFeature.SetField("city", sName)
                pFeature.SetField("latitude", dLat)
                pFeature.SetField("longitude", dLon)
                pFeature.SetField("population", dPop)
                pLayer_out.CreateFeature(pFeature)
                count += 1

    else:
        # Vector input dataset (Shapefile, GeoJSON, GeoPackage, etc.)
        pDataSource_in = ogr.Open(sFilename_in, 0)
        if pDataSource_in is None:
            raise RuntimeError(f"Could not open vector dataset: {sFilename_in}")

        pLayer_in = pDataSource_in.GetLayer(0)
        pSpatialRef_in = pLayer_in.GetSpatialRef()
        oTransform = None
        if pSpatialRef_in is not None and not pSpatialRef_in.IsSame(srs_target):
            oTransform = osr.CoordinateTransformation(pSpatialRef_in, srs_target)

        # Detect column names from layer definition
        layer_defn = pLayer_in.GetLayerDefn()
        field_names = [
            layer_defn.GetFieldDefn(i).GetName() for i in range(layer_defn.GetFieldCount())
        ]

        col_city_name = (
            sCol_city
            if sCol_city
            else next(
                (
                    c
                    for c in [
                        "NAME",
                        "name",
                        "CITY",
                        "city",
                        "NAMEASCII",
                        "nameascii",
                        "PLACE",
                        "place",
                    ]
                    if c in field_names
                ),
                None,
            )
        )
        col_pop_name = (
            sCol_pop
            if sCol_pop
            else next(
                (
                    c
                    for c in ["POP_MAX", "pop_max", "POPULATION", "population", "POP", "pop"]
                    if c in field_names
                ),
                None,
            )
        )

        pLayer_in.ResetReading()
        for pFeature_in in pLayer_in:
            pGeom = pFeature_in.GetGeometryRef()
            if pGeom is None:
                continue

            pGeom_target = pGeom.Clone()
            if oTransform is not None:
                pGeom_target.Transform(oTransform)

            dPop = (
                float(pFeature_in.GetField(col_pop_name))
                if col_pop_name and pFeature_in.GetField(col_pop_name) is not None
                else 0.0
            )
            if dPop < dPopulation_min:
                continue

            sName = (
                str(pFeature_in.GetField(col_city_name))
                if col_city_name and pFeature_in.GetField(col_city_name) is not None
                else ""
            )

            # Extract representative coordinate point
            if pGeom_target.GetGeometryType() == ogr.wkbPoint:
                dLon = pGeom_target.GetX()
                dLat = pGeom_target.GetY()
                pPoint = pGeom_target
            else:
                centroid = pGeom_target.Centroid()
                dLon = centroid.GetX()
                dLat = centroid.GetY()
                pPoint = centroid

            pFeature_out = ogr.Feature(pLayerDefn_out)
            pFeature_out.SetGeometry(pPoint)
            pFeature_out.SetField("city", sName)
            pFeature_out.SetField("latitude", dLat)
            pFeature_out.SetField("longitude", dLon)
            pFeature_out.SetField("population", dPop)
            pLayer_out.CreateFeature(pFeature_out)
            count += 1

        pDataSource_in = None

    pDataSource_out = None
    logger.info("Saved %d cities to %s", count, sFilename_out)
    return sFilename_out
