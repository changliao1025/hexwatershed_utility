"""
hexwatershed_utility
====================
A collection of tools for HexWatershed pre- and post-processing.
"""

__version__ = "0.1.0"
__author__ = "Chang Liao"

from . import core, postprocess, preprocess
from .core.exceptions import (
    DataInconsistencyError,
    HexWatershedError,
    MeshValidationError,
    SpatialReferenceError,
    TopologicalLoopError,
)

__all__ = [
    "core",
    "preprocess",
    "postprocess",
    "HexWatershedError",
    "MeshValidationError",
    "TopologicalLoopError",
    "DataInconsistencyError",
    "SpatialReferenceError",
]
