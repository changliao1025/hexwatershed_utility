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
    Rasterize lake vector polygons into a global/regional lake mask raster.

    Parameters
    ----------
    sFilename_in : str
        Path to input lake vector file (GeoJSON, Shapefile, etc.).
    sFilename_tif_out : str
        Path to output raster GeoTIFF file.
    dResolution_x_in : float, default=30/3600.0
        X-resolution in degrees (arc-seconds converted to degrees).
    dResolution_y_in : float, default=30/3600.0
        Y-resolution in degrees.
    iFlag_boundary_only_in : int, default=0
        If 1, rasterize boundary outline only. If 0, rasterize full polygon interior.

    Returns
    -------
    str
        Path to generated raster GeoTIFF.
    """
    if not os.path.exists(sFilename_in):
        raise FileNotFoundError(f"Input lake vector file not found: {sFilename_in}")

    sDir_out = os.path.dirname(sFilename_tif_out)
    if sDir_out and not os.path.exists(sDir_out):
        os.makedirs(sDir_out, exist_ok=True)

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
