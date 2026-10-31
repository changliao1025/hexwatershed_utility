import logging
import os

logger = logging.getLogger(__name__)


def create_city_mask(
    sFilename_in: str,
    sFilename_tif_out: str,
    dResolution_x_in: float = 30 / 3600.0,
    dResolution_y_in: float = 30 / 3600.0,
) -> str:
    """
    Rasterize city vector points into a global or regional raster mask.

    Parameters
    ----------
    sFilename_in : str
        Path to input city vector file (GeoJSON, Shapefile, etc.).
    sFilename_tif_out : str
        Path to output raster GeoTIFF file.
    dResolution_x_in : float, default=30/3600.0
        X-resolution in degrees.
    dResolution_y_in : float, default=30/3600.0
        Y-resolution in degrees.

    Returns
    -------
    str
        Path to output raster GeoTIFF file.
    """
    if not os.path.exists(sFilename_in):
        raise FileNotFoundError(f"Input city vector file not found: {sFilename_in}")

    sDir_out = os.path.dirname(sFilename_tif_out)
    if sDir_out and not os.path.exists(sDir_out):
        os.makedirs(sDir_out, exist_ok=True)

    from pyearth.toolbox.conversion.convert_vector_to_global_raster import (
        convert_vector_to_global_raster,
    )

    convert_vector_to_global_raster(
        sFilename_in, sFilename_tif_out, dResolution_x_in, dResolution_y_in
    )

    logger.info("Created city mask raster: %s", sFilename_tif_out)
    return sFilename_tif_out
