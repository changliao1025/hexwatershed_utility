import logging
import os

logger = logging.getLogger(__name__)


def create_lake_mask(
    sFilename_in: str,
    sFilename_tif_out: str,
    dResolution_x_in: float = 30 / 3600.0,
    dResolution_y_in: float = 30 / 3600.0,
    iFlag_boundary_only_in: int = 0,
) -> str:
    """
    Rasterize lake vector polygons into a global lake mask raster.

    Parameters
    ----------
    sFilename_in : str
        Input GeoJSON lake vector path.
    sFilename_tif_out : str
        Output raster GeoTIFF path.
    dResolution_x_in : float, default=30/3600.0
        X-resolution in degrees.
    dResolution_y_in : float, default=30/3600.0
        Y-resolution in degrees.
    iFlag_boundary_only_in : int, default=0
        Whether to rasterize boundary only (1) or entire polygon (0).

    Returns
    -------
    str
        Path to output raster file.
    """
    if not os.path.exists(sFilename_in):
        raise FileNotFoundError(f"Input lake vector not found: {sFilename_in}")

    from pyearth.toolbox.conversion.convert_vector_to_global_raster import (
        convert_vector_to_global_raster,
    )

    convert_vector_to_global_raster(
        sFilename_in,
        sFilename_tif_out,
        dResolution_x_in,
        dResolution_y_in,
        iFlag_boundary_only_in=iFlag_boundary_only_in,
    )

    logger.info("Created lake mask raster: %s", sFilename_tif_out)
    return sFilename_tif_out


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 2:
        create_lake_mask(sys.argv[1], sys.argv[2])
