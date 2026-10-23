import csv
import logging
import os
from osgeo import ogr, osr

logger = logging.getLogger(__name__)


def create_city_point(
    sFilename_csv_in: str,
    sFilename_geojson_out: str,
    dPopulation_min: float = 1.0e6,
) -> str:
    """
    Extract major cities exceeding a population threshold from a CSV file into GeoJSON points.

    Parameters
    ----------
    sFilename_csv_in : str
        Path to world cities CSV dataset.
    sFilename_geojson_out : str
        Path to output GeoJSON file.
    dPopulation_min : float, default=1.0e6
        Minimum population cutoff.

    Returns
    -------
    str
        Path to created GeoJSON file.
    """
    if not os.path.exists(sFilename_csv_in):
        raise FileNotFoundError(f"Input city CSV not found: {sFilename_csv_in}")

    pDriver = ogr.GetDriverByName("GeoJSON")
    if os.path.exists(sFilename_geojson_out):
        pDriver.DeleteDataSource(sFilename_geojson_out)

    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)

    pDataSource = pDriver.CreateDataSource(sFilename_geojson_out)
    pLayer = pDataSource.CreateLayer("city", srs, ogr.wkbPoint)

    pLayer.CreateField(ogr.FieldDefn("city", ogr.OFTString))
    pLayer.CreateField(ogr.FieldDefn("latitude", ogr.OFTReal))
    pLayer.CreateField(ogr.FieldDefn("longitude", ogr.OFTReal))
    pLayer.CreateField(ogr.FieldDefn("population", ogr.OFTReal))

    pLayerDefn = pLayer.GetLayerDefn()
    pFeature = ogr.Feature(pLayerDefn)

    count = 0
    with open(sFilename_csv_in, "r", encoding="utf-8", errors="ignore") as f:
        reader = csv.reader(f)
        header = next(reader, None)
        # Identify columns
        city_idx = 0
        lat_idx = 2
        lon_idx = 3
        pop_idx = 9

        for row in reader:
            if len(row) <= max(city_idx, lat_idx, lon_idx, pop_idx):
                continue
            sName = row[city_idx]
            try:
                dLatitude = float(row[lat_idx])
                dLongitude = float(row[lon_idx])
                pop_str = row[pop_idx].strip()
                if not pop_str:
                    continue
                dPopulation = float(pop_str)
            except ValueError:
                continue

            if dPopulation < dPopulation_min:
                continue

            pPoint = ogr.Geometry(ogr.wkbPoint)
            pPoint.AddPoint(dLongitude, dLatitude)
            pFeature.SetGeometry(pPoint)
            pFeature.SetField("city", sName)
            pFeature.SetField("latitude", dLatitude)
            pFeature.SetField("longitude", dLongitude)
            pFeature.SetField("population", dPopulation)
            pLayer.CreateFeature(pFeature)
            count += 1

    pDataSource = None
    logger.info("Saved %d cities to %s", count, sFilename_geojson_out)
    return sFilename_geojson_out


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 2:
        create_city_point(sys.argv[1], sys.argv[2])
