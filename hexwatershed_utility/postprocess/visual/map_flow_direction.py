import logging
import os

logger = logging.getLogger(__name__)


def map_flow_direction(
    sFilename_polyline_in: str,
    sFilename_output_in: str = "mpas_flow_direction_map.png",
    iFlag_global_in: int = 1,
) -> str:
    """
    Render 2D map of flow direction polylines.

    Parameters
    ----------
    sFilename_polyline_in : str
        Input GeoJSON polyline file of flow directions.
    sFilename_output_in : str, default='mpas_flow_direction_map.png'
        Output PNG image filename.
    iFlag_global_in : int, default=1
        Flag for global map projection.

    Returns
    -------
    str
        Path to output PNG.
    """
    if not os.path.exists(sFilename_polyline_in):
        raise FileNotFoundError(f"Polyline file not found: {sFilename_polyline_in}")

    from pyearthviz.map.vector.map_vector_polyline_file import map_vector_polyline_file

    map_vector_polyline_file(
        sFilename_polyline_in,
        iFlag_global_in=iFlag_global_in,
        sFilename_output_in=sFilename_output_in,
    )

    logger.info("Generated flow direction map: %s", sFilename_output_in)
    return sFilename_output_in


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        out = sys.argv[2] if len(sys.argv) > 2 else "mpas_flow_direction_map.png"
        map_flow_direction(sys.argv[1], out)
