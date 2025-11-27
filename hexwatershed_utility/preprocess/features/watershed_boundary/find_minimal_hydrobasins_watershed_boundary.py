import os
import glob
from contextlib import contextmanager
from typing import Optional, List, Union, Dict, Any
from osgeo import gdal, ogr
gdal.UseExceptions()
from pyearth.system.define_global_variables import *
from pyearth.gis.gdal.read.vector.gdal_get_vector_extent import gdal_get_vector_extent
from pyearth.gis.location.get_hydrosheds_continent_from_extent import get_hydrosheds_continent_from_extent
from hexwatershed_utility.preprocess.features.watershed_boundary.shapefile_repair_utility import open_shapefile_with_fallback, diagnose_and_fix_shapefile



def _find_containing_watershed(pGeometry_network: ogr.Geometry,
                              sFilename_watershed: str) -> Optional[str]:
    """
    Find watershed polygon that contains the given geometry.

    Args:
        pGeometry_network: OGR Geometry to search for
        sFilename_watershed: Path to watershed shapefile

    Returns:
        WKT string of containing watershed geometry, or None if not found
    """
    try:
        # Use the shapefile repair utility to handle corrupted files
        pDataset_watershed = ogr.Open(sFilename_watershed, 0)
        if pDataset_watershed is None:
            print(f"Could not open watershed file even with repair attempts: {sFilename_watershed}")
            return None

        pLayer_watershed = pDataset_watershed.GetLayer(0)
        if pLayer_watershed is None:
            print(f"Could not get layer from watershed file: {sFilename_watershed}")
            pDataset_watershed = None
            return None

        # Set spatial filter to improve performance
        aEnvelope = pGeometry_network.GetEnvelope()
        try:
            aExtent = pLayer_watershed.GetExtent()
            minx, maxx, miny, maxy = aExtent
        except Exception as e:
            print(f"Could not get layer extent, using geometry envelope: {e}")

        # Try to set spatial filter - this can fail with "Inconsistent shape count for bin"

        pLayer_watershed.ResetReading()
        #nFeature_count = pLayer_watershed.GetFeatureCount()
        #for i in range(nFeature_count):
        for pFeature_watershed in pLayer_watershed:
            pGeometry_watershed = pFeature_watershed.GetGeometryRef()
            if pGeometry_watershed is None:
                continue

            # Check containment first (more restrictive)
            try:
                if pGeometry_watershed.Contains(pGeometry_network):
                    result = pGeometry_watershed.ExportToWkt()
                    pDataset_watershed = None
                    return result
            except Exception as e:
                print(f"Error checking containment: {e}")
                continue

            # Check intersection (less restrictive - causes skip to next level)
            try:
                if pGeometry_watershed.Intersects(pGeometry_network):
                    print(f"River network intersects with watershed, skipping to next level")
                    pDataset_watershed = None
                    return "INTERSECTS"  # Special flag to indicate intersection
            except Exception as e:
                print(f"Error checking intersection: {e}")
                continue

        pDataset_watershed = None

    except Exception as e:
        error_msg = str(e)
        print(f"Error processing watershed file {sFilename_watershed}: {error_msg}")

        # Check if this is the specific "Inconsistent shape count for bin" error
        if "Inconsistent shape count for bin" in error_msg:
            print("Detected shapefile corruption. Attempting diagnosis and repair...")
            success, message = diagnose_and_fix_shapefile(sFilename_watershed)
            if success:
                print(f"Repair successful: {message}")
                print("Retrying with repaired shapefile...")
                # Retry the operation once after repair
                try:
                    return _find_containing_watershed(pGeometry_network, sFilename_watershed)
                except Exception as retry_e:
                    print(f"Retry failed even after repair: {retry_e}")
            else:
                print(f"Repair failed: {message}")

    return None

def _get_watershed_shapefiles(sFolder: str, iLevel: int) -> List[str]:
    """Get watershed shapefile paths for a given level."""
    sLevel = f"lev{iLevel:02d}"
    sFilename_reg = f'hybas_lake_*{sLevel}*_*.shp'
    return glob.glob(os.path.join(sFolder, sFilename_reg))

def _find_watershed_for_geometry(pGeometry_network: ogr.Geometry, sFolder: str,
                               iLevel_start: int, iLevel_end: int) -> Optional[str]:
    """Find the minimal watershed boundary for a given geometry using your naming style."""
    for iLevel in range(iLevel_start, iLevel_end, -1):
        aFilename = _get_watershed_shapefiles(sFolder, iLevel)
        if not aFilename:
            print(f"No shapefiles found for level lev{iLevel:02d} in {sFolder}")
            continue

        sFilename_full = aFilename[0]  # Use the first matching shapefile
        sResult = _find_containing_watershed(pGeometry_network, sFilename_full)

        if sResult == "INTERSECTS":
            # Skip to next level (larger watersheds)
            continue
        elif sResult is not None:
            # Found containing watershed
            print(f"Found containing watershed at level lev{iLevel:02d}")
            return sResult

    print("No containing watershed boundary found")
    return None

def find_minimal_hydrobasins_watershed_boundary(sFilename_river_network_in: str,
                                               sFolder_watershed_boundary_in: str,
                                               iFlag_nested_in: bool = False) -> Optional[Union[str, List[str]]]:
    """
    Find the minimal watershed boundary that completely contains a river network.

    Searches from low level (small watersheds) to high level (large watersheds)
    to find the smallest boundary that encompasses the entire river network.

    Args:
        sFilename_river_network_in (str): Path to river network GeoJSON file
        sFolder_watershed_boundary_in (str): Path to folder containing watershed boundary shapefiles
                                           organized by levels (e.g., level_01, level_02, etc.)

    Returns:
        dict: Dictionary containing:
            - 'level': The watershed level (e.g., 'level_03')
            - 'filename': Path to the shapefile containing the boundary
            - 'feature_id': ID/index of the feature that contains the network
            - 'geometry': The boundary geometry as OGR geometry object
            - 'bounds': Bounding box of the containing watershed (minx, miny, maxx, maxy)
            - 'area': Area of the containing watershed
    """



    # Main function logic
    print(f"Reading river network from: {sFilename_river_network_in}")

    # Get river network bounds and geometry
    aNetwork_bounds, pGeometry_network = gdal_get_vector_extent(sFilename_river_network_in, iFlag_return_union_geometry=1)
    if aNetwork_bounds is None:
        print("Failed to read river network")
        return None

    #get containing folder of the river network file
    sFolder_river_network = os.path.dirname(sFilename_river_network_in)

    print(f"River network bounds: {aNetwork_bounds}")

    # use regex to find the folder what has the name of the region such as 'na' for north america

    sRegion = get_hydrosheds_continent_from_extent(aNetwork_bounds)
    sRegax = '*_' + sRegion + '_*'

    aFolder = glob.glob(os.path.join(sFolder_watershed_boundary_in, sRegax))
    if not aFolder:
        print(f"No watershed folders found for region '{sRegion}' in {sFolder_watershed_boundary_in}")
        return None

    for sFolder in aFolder:
        if not os.path.isdir(sFolder):
            print(f"Skipping non-directory: {sFolder}")
            continue
        else:
            print(f"Searching in watershed folder: {sFolder}")
            break

    iLevel_start = 10
    iLevel_end = 1

    # Process based on nested flag
    if iFlag_nested_in:
        # Process each feature separately
        aWatershed_boundaries = list()

        try:
            pDataset_network = ogr.Open(sFilename_river_network_in, 0)

            if pDataset_network is None:
                print(f"Could not open river network file: {sFilename_river_network_in}")
                return None
            pLayer_network = pDataset_network.GetLayer(0)
            if pLayer_network is None:
                print(f"Could not get layer from river network file: {sFilename_river_network_in}")
                return None
            nSubbasin = pLayer_network.GetFeatureCount()
            if nSubbasin == 0:
                print(f"No features found in river network file: {sFilename_river_network_in}")
                return None
            pLayer_network.ResetReading()
            for pFeature_network in pLayer_network:
            #for i in range(0, 10, 1):
                #pFeature_network = pLayer_network.GetFeature(i)
                pGeometry_network = pFeature_network.GetGeometryRef()
                if pGeometry_network is None:
                    print("Feature has no valid geometry, skipping")
                    continue
                # Find watershed for this feature using optimized helper function
                sWkt = _find_watershed_for_geometry(pGeometry_network, sFolder, iLevel_start, iLevel_end)
                if sWkt and sWkt != "INTERSECTS":
                    aWatershed_boundaries.append(sWkt)
                    print(f"Watershed boundary for subbasin {pFeature_network.GetFID()} found")

        except Exception as e:
            print(f"Error processing nested features: {e}")
            return None

        # Return individual watershed boundaries
        if len(aWatershed_boundaries) == 0:
            print("No watershed boundaries found for any subbasin")
            return None
        else:
            print(f"Found {len(aWatershed_boundaries)} individual watershed boundaries")
            return aWatershed_boundaries
    else:
        # Process union geometry of all features
        sWkt = _find_watershed_for_geometry(pGeometry_network, sFolder, iLevel_start, iLevel_end)
        if sWkt and sWkt != "INTERSECTS":
            return sWkt

        print("No containing watershed boundary found for", pGeometry_network.ExportToWkt())
        return None
