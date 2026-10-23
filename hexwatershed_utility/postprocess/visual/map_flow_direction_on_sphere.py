import logging
import os

logger = logging.getLogger(__name__)


def map_flow_direction_on_sphere(sFilename_polyline_in: str):
    """
    Render 3D animated sphere of flow direction polylines using GeoVista.

    Parameters
    ----------
    sFilename_polyline_in : str
        Input GeoJSON polyline file of flow directions.
    """
    if not os.path.exists(sFilename_polyline_in):
        raise FileNotFoundError(f"Polyline file not found: {sFilename_polyline_in}")

    from pyearthviz3d.geovista.animate_polyline_file_on_sphere import (
        animate_polyline_file_on_sphere,
    )

    return animate_polyline_file_on_sphere(sFilename_polyline_in)


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        map_flow_direction_on_sphere(sys.argv[1])
