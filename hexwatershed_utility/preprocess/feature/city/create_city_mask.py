import logging
import os

logger = logging.getLogger(__name__)


def create_city_mask(
    sFilename_geojson_in: str,
    sFilename_tif_out: str,
    dResolution_x_in: float = 30 / 3600.0,
    dResolution_y_in: float = 30 / 3600.0,
) -> str:
    """
    Rasterize city vector points into a global raster mask.

    Parameters
    ----------
    sFilename_geojson_in : str
        Input GeoJSON point file.
    sFilename_tif_out : str
        Output raster GeoTIFF file.
    dResolution_x_in : float, default=30/3600.0
        X-resolution in degrees.
    dResolution_y_in : float, default=30/3600.0
        Y-resolution in degrees.

    Returns
    -------
    str
        Path to output raster.
    """
    if not os.path.exists(sFilename_geojson_in):
        raise FileNotFoundError(f"Input city GeoJSON not found: {sFilename_geojson_in}")

    from pyearth.toolbox.conversion.convert_vector_to_global_raster import (
        convert_vector_to_global_raster,
    )

    convert_vector_to_global_raster(
        sFilename_geojson_in, sFilename_tif_out, dResolution_x_in, dResolution_y_in
    )

    logger.info("Created city mask raster: %s", sFilename_tif_out)
    return sFilename_tif_out


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 2:
        create_city_mask(sys.argv[1], sys.argv[2])
