import json
import logging
import requests

logger = logging.getLogger(__name__)

URL_CONUS = "https://services1.arcgis.com/cRvLdSPAsRupRo7I/arcgis/rest/services/US_State_Boundaries_(CONUS)/FeatureServer/0/query"


def download_conus_boundary(sFilename_geojson_out: str = "US_State_Boundaries_CONUS.geojson") -> str:
    """
    Download CONUS state boundaries from ArcGIS REST service as GeoJSON.

    Parameters
    ----------
    sFilename_geojson_out : str, default='US_State_Boundaries_CONUS.geojson'
        Output GeoJSON path.

    Returns
    -------
    str
        Path to the saved GeoJSON file.
    """
    params = {
        "where": "1=1",
        "outFields": "*",
        "f": "geojson",
    }

    response = requests.get(URL_CONUS, params=params, timeout=60)
    if response.status_code == 200:
        with open(sFilename_geojson_out, "w") as f:
            json.dump(response.json(), f)
        logger.info("Data downloaded and saved to %s", sFilename_geojson_out)
        return sFilename_geojson_out
    else:
        raise RuntimeError(
            f"Failed to download CONUS boundary. HTTP status code: {response.status_code}"
        )


if __name__ == "__main__":
    download_conus_boundary()
