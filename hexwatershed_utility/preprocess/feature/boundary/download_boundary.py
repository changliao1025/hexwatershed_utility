import json
import logging
import os
from typing import Optional

import requests

logger = logging.getLogger(__name__)

URL_CONUS_DEFAULT = (
    "https://services1.arcgis.com/cRvLdSPAsRupRo7I/arcgis/rest/services/"
    "US_State_Boundaries_(CONUS)/FeatureServer/0/query"
)


def download_boundary(
    sUrl: Optional[str] = None,
    sFilename_geojson_out: str = "boundary.geojson",
    sWhere: str = "1=1",
    sOutFields: str = "*",
    timeout: int = 60,
) -> str:
    """
    Download administrative or regional boundaries from an ArcGIS REST service as GeoJSON.

    Parameters
    ----------
    sUrl : Optional[str], default=None
        ArcGIS REST query endpoint URL. If None, defaults to CONUS state boundaries service.
    sFilename_geojson_out : str, default='boundary.geojson'
        Path to save output GeoJSON file.
    sWhere : str, default='1=1'
        SQL where clause for the query.
    sOutFields : str, default='*'
        Comma-separated list of fields to include.
    timeout : int, default=60
        Request timeout in seconds.

    Returns
    -------
    str
        Path to downloaded GeoJSON file.
    """
    url = sUrl or URL_CONUS_DEFAULT
    params = {
        "where": sWhere,
        "outFields": sOutFields,
        "f": "geojson",
    }

    sDir = os.path.dirname(sFilename_geojson_out)
    if sDir and not os.path.exists(sDir):
        os.makedirs(sDir, exist_ok=True)

    logger.info("Querying boundary from: %s", url)
    response = requests.get(url, params=params, timeout=timeout)
    if response.status_code == 200:
        data = response.json()
        with open(sFilename_geojson_out, "w", encoding="utf-8") as f:
            json.dump(data, f)
        logger.info("Saved boundary GeoJSON to %s", sFilename_geojson_out)
        return sFilename_geojson_out
    else:
        raise RuntimeError(f"Failed to download boundary. HTTP status code: {response.status_code}")


def download_conus_boundary(
    sFilename_geojson_out: str = "US_State_Boundaries_CONUS.geojson",
) -> str:
    """Convenience shortcut to download CONUS state boundaries."""
    return download_boundary(sUrl=URL_CONUS_DEFAULT, sFilename_geojson_out=sFilename_geojson_out)
