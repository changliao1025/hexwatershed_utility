"""
Shapefile Repair Utility

This module provides functions to handle and repair corrupted shapefiles that cause
"Inconsistent shape count for bin" errors in GDAL/OGR operations.
"""

import os
import shutil
from typing import Optional, Tuple
from osgeo import gdal, ogr

def diagnose_shapefile_corruption(sFilename_shapefile: str) -> dict:
    """
    Diagnose shapefile corruption issues.

    Args:
        sFilename_shapefile: Path to shapefile (.shp)

    Returns:
        Dictionary with diagnostic information
    """
    diagnostic_info = {
        'shp_exists': False,
        'shx_exists': False,
        'dbf_exists': False,
        'prj_exists': False,
        'shx_consistent': False,
        'can_open': False,
        'feature_count_error': False,
        'geometry_errors': 0,
        'recommended_action': 'unknown'
    }

    base_path = os.path.splitext(sFilename_shapefile)[0]

    # Check file existence
    diagnostic_info['shp_exists'] = os.path.exists(sFilename_shapefile)
    diagnostic_info['shx_exists'] = os.path.exists(base_path + '.shx')
    diagnostic_info['dbf_exists'] = os.path.exists(base_path + '.dbf')
    diagnostic_info['prj_exists'] = os.path.exists(base_path + '.prj')

    if not diagnostic_info['shp_exists']:
        diagnostic_info['recommended_action'] = 'file_missing'
        return diagnostic_info

    # Check if we can open the shapefile
    try:
        pDataset = ogr.Open(sFilename_shapefile, 0)
        if pDataset is not None:
            diagnostic_info['can_open'] = True
            pLayer = pDataset.GetLayer(0)

            # Try to get feature count - this often triggers the error
            try:
                nFeatures = pLayer.GetFeatureCount()
                diagnostic_info['feature_count'] = nFeatures
            except Exception as e:
                diagnostic_info['feature_count_error'] = True
                diagnostic_info['feature_count_error_msg'] = str(e)

            # Check geometry validity for first few features
            pLayer.ResetReading()
            nChecked = 0
            nGeomErrors = 0
            for pFeature in pLayer:
                if nChecked >= 10:  # Check first 10 features
                    break
                pGeometry = pFeature.GetGeometryRef()
                if pGeometry is None or not pGeometry.IsValid():
                    nGeomErrors += 1
                nChecked += 1

            diagnostic_info['geometry_errors'] = nGeomErrors
            diagnostic_info['features_checked'] = nChecked
            pDataset = None

    except Exception as e:
        diagnostic_info['open_error'] = str(e)

    # Determine recommended action
    if diagnostic_info['feature_count_error']:
        if not diagnostic_info['shx_exists']:
            diagnostic_info['recommended_action'] = 'rebuild_shx'
        else:
            diagnostic_info['recommended_action'] = 'repair_shx'
    elif diagnostic_info['geometry_errors'] > 0:
        diagnostic_info['recommended_action'] = 'geometry_repair'
    elif diagnostic_info['can_open']:
        diagnostic_info['recommended_action'] = 'no_action_needed'
    else:
        diagnostic_info['recommended_action'] = 'file_corrupted'

    return diagnostic_info

def repair_shapefile_spatial_index(sFilename_shapefile: str, sBackup_suffix: str = '_backup') -> bool:
    """
    Repair shapefile spatial index (.shx file) by regenerating it.

    Args:
        sFilename_shapefile: Path to shapefile
        sBackup_suffix: Suffix for backup files

    Returns:
        True if repair succeeded, False otherwise
    """
    try:
        base_path = os.path.splitext(sFilename_shapefile)[0]
        shx_file = base_path + '.shx'

        # Create backup of original .shx if it exists
        if os.path.exists(shx_file):
            backup_shx = base_path + sBackup_suffix + '.shx'
            shutil.copy2(shx_file, backup_shx)
            print(f"Backed up {shx_file} to {backup_shx}")

        # Remove corrupted .shx file
        if os.path.exists(shx_file):
            os.remove(shx_file)
            print(f"Removed corrupted {shx_file}")

        # Set GDAL options to rebuild spatial index
        gdal.SetConfigOption('SHAPE_RESTORE_SHX', 'YES')

        # Open shapefile - this should trigger .shx rebuild
        pDataset = ogr.Open(sFilename_shapefile, 0)
        if pDataset is None:
            print("Failed to open shapefile for repair")
            return False

        pLayer = pDataset.GetLayer(0)
        if pLayer is None:
            print("Failed to access layer for repair")
            return False

        # Force reading through all features to rebuild index
        pLayer.ResetReading()
        nFeatures = 0
        for pFeature in pLayer:
            nFeatures += 1

        pDataset = None
        print(f"Successfully rebuilt spatial index. Found {nFeatures} features")

        # Verify repair worked
        pDataset_test = ogr.Open(sFilename_shapefile, 0)
        if pDataset_test is not None:
            pLayer_test = pDataset_test.GetLayer(0)
            try:
                nCount = pLayer_test.GetFeatureCount()
                pDataset_test = None
                print(f"Repair verification: {nCount} features accessible")
                return True
            except Exception as e:
                print(f"Repair verification failed: {e}")
                pDataset_test = None
                return False

        return False

    except Exception as e:
        print(f"Shapefile repair failed: {e}")
        return False
    finally:
        gdal.SetConfigOption('SHAPE_RESTORE_SHX', None)

def open_shapefile_with_fallback(sFilename_shapefile: str) -> Optional[ogr.DataSource]:
    """
    Open shapefile with multiple fallback strategies for corrupted files.

    Args:
        sFilename_shapefile: Path to shapefile

    Returns:
        OGR DataSource or None if all methods failed
    """
    print(f"Attempting to open: {sFilename_shapefile}")

    # Strategy 1: Normal opening
    try:
        pDataset = ogr.Open(sFilename_shapefile, 0)
        if pDataset is not None:
            pLayer = pDataset.GetLayer(0)
            if pLayer is not None:
                # Quick test
                pLayer.ResetReading()
                nCount = pLayer.GetFeatureCount()
                print(f"Successfully opened normally ({nCount} features)")
                return pDataset
    except Exception as e:
        print(f"Normal opening failed: {e}")
        if pDataset:
            pDataset = None

    # Strategy 2: Force spatial index rebuild
    try:
        gdal.SetConfigOption('SHAPE_RESTORE_SHX', 'YES')
        pDataset = ogr.Open(sFilename_shapefile, 0)
        if pDataset is not None:
            pLayer = pDataset.GetLayer(0)
            if pLayer is not None:
                pLayer.ResetReading()
                print("Successfully opened with spatial index rebuild")
                return pDataset
    except Exception as e:
        print(f"Spatial index rebuild failed: {e}")
        if pDataset:
            pDataset = None
    finally:
        gdal.SetConfigOption('SHAPE_RESTORE_SHX', None)

    # Strategy 3: Disable spatial index
    try:
        gdal.SetConfigOption('SHAPE_RESTORE_SHX', 'NO')
        pDataset = ogr.Open(sFilename_shapefile, 0)
        if pDataset is not None:
            pLayer = pDataset.GetLayer(0)
            if pLayer is not None:
                pLayer.ResetReading()
                print("Successfully opened without spatial index")
                return pDataset
    except Exception as e:
        print(f"No spatial index opening failed: {e}")
        if pDataset:
            pDataset = None
    finally:
        gdal.SetConfigOption('SHAPE_RESTORE_SHX', None)

    # Strategy 4: Manual repair attempt
    print("Attempting manual repair of spatial index...")
    if repair_shapefile_spatial_index(sFilename_shapefile):
        try:
            pDataset = ogr.Open(sFilename_shapefile, 0)
            if pDataset is not None:
                pLayer = pDataset.GetLayer(0)
                if pLayer is not None:
                    print("Successfully opened after manual repair")
                    return pDataset
        except Exception as e:
            print(f"Opening after repair failed: {e}")

    print("All opening strategies failed")
    return None

def get_safe_feature_iterator(pLayer: ogr.Layer, iMax_features: Optional[int] = None):
    """
    Get a safe iterator for features that handles corrupted geometries.

    Args:
        pLayer: OGR Layer
        iMax_features: Maximum number of features to iterate (None for all)

    Yields:
        Valid OGR Features
    """
    pLayer.ResetReading()
    nProcessed = 0

    while True:
        if iMax_features is not None and nProcessed >= iMax_features:
            break

        try:
            pFeature = pLayer.GetNextFeature()
            if pFeature is None:
                break

            # Check if feature has valid geometry
            pGeometry = pFeature.GetGeometryRef()
            if pGeometry is not None and pGeometry.IsValid():
                yield pFeature
                nProcessed += 1
            else:
                print(f"Skipping feature {pFeature.GetFID()} with invalid geometry")

        except Exception as e:
            print(f"Error reading feature: {e}")
            break

# Usage example and diagnostic function
def diagnose_and_fix_shapefile(sFilename_shapefile: str) -> Tuple[bool, str]:
    """
    Comprehensive diagnosis and repair of shapefile corruption.

    Args:
        sFilename_shapefile: Path to shapefile

    Returns:
        Tuple of (success, message)
    """
    print(f"\n=== Diagnosing shapefile: {sFilename_shapefile} ===")

    # Run diagnosis
    diagnostic = diagnose_shapefile_corruption(sFilename_shapefile)

    print("Diagnostic Results:")
    for key, value in diagnostic.items():
        if not key.endswith('_msg'):
            print(f"  {key}: {value}")

    # Attempt repair based on diagnosis
    if diagnostic['recommended_action'] == 'no_action_needed':
        return True, "Shapefile is healthy"
    elif diagnostic['recommended_action'] == 'file_missing':
        return False, "Shapefile does not exist"
    elif diagnostic['recommended_action'] in ['rebuild_shx', 'repair_shx']:
        print("\nAttempting spatial index repair...")
        if repair_shapefile_spatial_index(sFilename_shapefile):
            return True, "Spatial index successfully repaired"
        else:
            return False, "Spatial index repair failed"
    elif diagnostic['recommended_action'] == 'file_corrupted':
        return False, "Shapefile is severely corrupted and cannot be repaired"
    else:
        return False, f"Unknown issue: {diagnostic['recommended_action']}"