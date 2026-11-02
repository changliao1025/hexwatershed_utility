import logging
import os

from pyearth.toolbox.analysis.extract.clip_vector_by_polygon_file import clip_vector_by_polygon_file

logger = logging.getLogger(__name__)


def crop_river_network_by_boundary(
    sFilename_river: str,
    sFilename_boundary: str,
    sFilename_out: str,
    iFlag_endorheic: int = 0,
) -> None:
    """
    Crop a river network vector file by a watershed boundary polygon.

    Parameters
    ----------
    sFilename_river : str
        Path to the input river network vector dataset (e.g., .shp, .geojson, .gpkg).
    sFilename_boundary : str
        Path to the clipping boundary polygon vector dataset.
    sFilename_out : str
        Path to the output clipped river network vector dataset.
    iFlag_endorheic : int, default=0
        If 1, skip merging for endorheic basins. Default is 0.
    """
    if not os.path.exists(sFilename_river):
        raise FileNotFoundError(f"Input river network file not found: {sFilename_river}")
    if not os.path.exists(sFilename_boundary):
        raise FileNotFoundError(f"Clipping boundary file not found: {sFilename_boundary}")

    sFolder_out = os.path.dirname(sFilename_out)
    if sFolder_out:
        os.makedirs(sFolder_out, exist_ok=True)

    logger.info(
        "Cropping river network %s using boundary %s...", sFilename_river, sFilename_boundary
    )

    clip_vector_by_polygon_file(
        sFilename_vector_in=sFilename_river,
        sFilename_polygon_in=sFilename_boundary,
        sFilename_vector_out=sFilename_out,
        iFlag_endorheic=iFlag_endorheic,
    )

    logger.info("River network cropped and saved to: %s", sFilename_out)


__all__ = ["crop_river_network_by_boundary"]
