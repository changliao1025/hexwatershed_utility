"""
Create MOSART bifurcation NetCDF file with downstream dimension.
Supports both structured (2D->3D) and unstructured (1D->2D) mesh formats.
Adds bifurcation_ratio and ibt_demand variables.

Converted from MATLAB to Python.
"""

import numpy as np
import pandas as pd
import netCDF4 as nc
import shutil
from pathlib import Path
from typing import Tuple, Optional


def create_bifurc_netcdf(
    input_file: str,
    output_file: str,
    csv_file: Optional[str] = None,
    max_downstream: int = 2
) -> None:
    """
    Create MOSART bifurcation NetCDF file with downstream dimension.

    Parameters
    ----------
    input_file : str
        Path to input MOSART NetCDF file
    output_file : str
        Path to output NetCDF file with bifurcation support
    csv_file : str, optional
        Path to CSV file containing bifurcation/IBT data
    max_downstream : int, default=2
        Maximum number of downstream connections per cell
    """
    print(f'Reading original NetCDF file: {input_file}')

    # Open input file to detect mesh type
    with nc.Dataset(input_file, 'r') as ds:
        # Detect mesh type by checking cell ID dimensions
        ID_var = ds.variables['ID']
        ndims = len(ID_var.shape)
        is_structured = (ndims == 2)

        # Get the lat/lon for each cell
        lat = ds.variables['latixy'][:]
        lon = ds.variables['longxy'][:]
        ID = ds.variables['ID'][:]

        if is_structured:
            # Structured mesh: dnID(lon(row), lat(col))
            dnID_2d = ds.variables['dnID'][:]
            nrow, ncol = dnID_2d.shape
            ncells = nrow * ncol
            print(f'Structured mesh detected: dnID dimensions {nrow} x {ncol} (total {ncells} cells)')
        else:
            # Unstructured mesh: dnID(cells)
            dnID_1d = ds.variables['dnID'][:]
            ncells = len(dnID_1d)
            nrow = ncells  # For unstructured, treat as 1D array
            ncol = 1       # Dummy dimension
            print(f'Unstructured mesh detected: dnID dimensions {ncells} cells')

    # Create arrays with downstream dimension
    if is_structured:
        # Structured: Create 3D arrays
        dnID_2d_expand_to_3d = -9999 * np.ones((nrow, ncol, max_downstream), dtype=np.int32)
        bifurc_ratio_3d = np.zeros((nrow, ncol, max_downstream), dtype=np.float64)
        ibt_demand_3d = np.zeros((nrow, ncol, max_downstream), dtype=np.float64)

        # Copy original 2D data to first layer (primary downstream)
        dnID_2d_expand_to_3d[:, :, 0] = dnID_2d
    else:
        # Unstructured: Create 2D arrays
        dnID_1d_expand_to_2d = -9999 * np.ones((ncells, max_downstream), dtype=np.int32)
        bifurc_ratio_2d = np.zeros((ncells, max_downstream), dtype=np.float64)
        ibt_demand_2d = np.zeros((ncells, max_downstream), dtype=np.float64)

        # Copy original 1D data to first column (primary downstream)
        dnID_1d_expand_to_2d[:, 0] = dnID_1d

    # Copy original file to new file
    print(f'Copying original file to: {output_file}')
    shutil.copyfile(input_file, output_file)

    # Open the copied file for modification
    with nc.Dataset(output_file, 'r+') as ds:
        # Create downstream dimension
        ds.createDimension('downstream', max_downstream)

        # Get original dnID variable info
        dnID_var = ds.variables['dnID']
        dnID_dtype = dnID_var.dtype

        if is_structured:
            # Structured mesh: Get existing dimension names
            lon_dim = 'lon'
            lat_dim = 'lat'

            # Rename original dnID variable to keep it
            ds.renameVariable('dnID', 'dnID_original')

            # Create new 3D variables: dnID(lon, lat, downstream)
            dnID_new = ds.createVariable('dnID', dnID_dtype, (lon_dim, lat_dim, 'downstream'))
            bifurc_ratio_var = ds.createVariable('bifurc_ratio', 'f8', (lon_dim, lat_dim, 'downstream'))
            ibt_demand_var = ds.createVariable('ibt_demand', 'f8', (lon_dim, lat_dim, 'downstream'))

            # Add variable attributes
            bifurc_ratio_var.long_name = 'Bifurcation split ratios'
            bifurc_ratio_var.units = 'dimensionless'
            bifurc_ratio_var.valid_range = np.array([0.0, 1.0])

            ibt_demand_var.long_name = 'Inter-basin transfer demand'
            ibt_demand_var.units = 'm3/s'
            ibt_demand_var.valid_min = 0.0

        else:
            # Unstructured mesh: Get existing dimension name
            cells_dim = 'gridcell'  # Standard MPAS name

            # Rename original dnID variable to keep it
            ds.renameVariable('dnID', 'dnID_original')

            # Create new 2D variables: dnID(cells, downstream)
            dnID_new = ds.createVariable('dnID', dnID_dtype, (cells_dim, 'downstream'))
            bifurc_ratio_var = ds.createVariable('bifurc_ratio', 'f8', (cells_dim, 'downstream'))
            ibt_demand_var = ds.createVariable('ibt_demand', 'f8', (cells_dim, 'downstream'))

            # Add variable attributes
            bifurc_ratio_var.long_name = 'Bifurcation split ratios'
            bifurc_ratio_var.units = 'dimensionless'
            bifurc_ratio_var.valid_range = np.array([0.0, 1.0])

            ibt_demand_var.long_name = 'Inter-basin transfer demand'
            ibt_demand_var.units = 'm3/s'
            ibt_demand_var.valid_min = 0.0

        # Write the data
        if is_structured:
            print('Writing 3D structured mesh data...')
            dnID_new[:] = dnID_2d_expand_to_3d
            bifurc_ratio_var[:] = bifurc_ratio_3d
            ibt_demand_var[:] = ibt_demand_3d
            print('Structured format: dnID(lon,lat,downstream), bifurc_ratio(lon,lat,downstream), ibt_demand(lon,lat,downstream)')
        else:
            print('Writing 2D unstructured mesh data...')
            dnID_new[:] = dnID_1d_expand_to_2d
            bifurc_ratio_var[:] = bifurc_ratio_2d
            ibt_demand_var[:] = ibt_demand_2d
            print('Unstructured format: dnID(cells,downstream), bifurc_ratio(cells,downstream), ibt_demand(cells,downstream)')

    print('Successfully created bifurcation NetCDF file!')
    print('Note: Original dnID preserved as "dnID_original"')

    # Process CSV file if provided
    if csv_file is not None:
        print(f'\nProcessing CSV file: {csv_file}')
        bif = pd.read_csv(csv_file)
        print(f'Found {len(bif)} bifurcation/IBT entries')

        # Initialize arrays for cell IDs to be filled in CSV
        cell_ids = np.zeros(len(bif), dtype=np.int32)
        dncell1_ids = np.zeros(len(bif), dtype=np.int32)
        dncell2_ids = np.zeros(len(bif), dtype=np.int32)

        # Process each row in the CSV
        for i in range(len(bif)):
            print(f'\nProcessing row {i}: {bif.iloc[i]["note"]}')

            # Extract data from CSV row
            split_lat = bif.iloc[i]['cell_lat']
            split_lon = bif.iloc[i]['cell_long']
            is_ibt = bif.iloc[i]['ibt']

            # Convert longitude from 0-360 to -180-180 format if needed
            if split_lon > 180:
                split_lon = split_lon - 360

            # Find closest cell for splitting point
            split_cell_idx, split_cell_id = find_closest_cell(
                lat, lon, ID, split_lat, split_lon, is_structured
            )
            if split_cell_idx == -1:
                raise ValueError(f'Could not find cell for splitting point at lat={split_lat:.4f}, lon={split_lon:.4f}')

            print(f'  Split point: lat={split_lat:.4f}, lon={split_lon:.4f} -> cell_id={split_cell_id}')
            cell_ids[i] = split_cell_id

            # Handle primary downstream cell
            if bif.iloc[i]['dncell1_lat'] != -9999 and bif.iloc[i]['dncell1_long'] != -9999:
                # Custom primary downstream specified
                dn1_lat = bif.iloc[i]['dncell1_lat']
                dn1_lon = bif.iloc[i]['dncell1_long']

                dn1_cell_idx, dn1_cell_id = find_closest_cell(
                    lat, lon, ID, dn1_lat, dn1_lon, is_structured
                )
                if dn1_cell_idx == -1:
                    raise ValueError(f'Could not find primary downstream cell at lat={dn1_lat:.4f}, lon={dn1_lon:.4f}')

                print(f'  Primary downstream: lat={dn1_lat:.4f}, lon={dn1_lon:.4f} -> cell_id={dn1_cell_id}')
                dncell1_ids[i] = dn1_cell_id

                # Update primary downstream in arrays
                if is_structured:
                    split_row, split_col = get_structured_indices(split_cell_idx, nrow, ncol)
                    dnID_2d_expand_to_3d[split_row, split_col, 0] = dn1_cell_id
                else:
                    dnID_1d_expand_to_2d[split_cell_idx, 0] = dn1_cell_id
            else:
                # Use existing primary downstream
                if is_structured:
                    split_row, split_col = get_structured_indices(split_cell_idx, nrow, ncol)
                    existing_dnID = dnID_2d_expand_to_3d[split_row, split_col, 0]
                else:
                    existing_dnID = dnID_1d_expand_to_2d[split_cell_idx, 0]

                dncell1_ids[i] = existing_dnID
                print(f'  Primary downstream: using existing cell_id={existing_dnID}')

            # Handle secondary downstream cell
            dn2_lat = bif.iloc[i]['dncell2_lat']
            dn2_lon = bif.iloc[i]['dncell2_long']

            dn2_cell_idx, dn2_cell_id = find_closest_cell(
                lat, lon, ID, dn2_lat, dn2_lon, is_structured
            )
            if dn2_cell_idx == -1:
                raise ValueError(f'Could not find secondary downstream cell at lat={dn2_lat:.4f}, lon={dn2_lon:.4f}')

            print(f'  Secondary downstream: lat={dn2_lat:.4f}, lon={dn2_lon:.4f} -> cell_id={dn2_cell_id}')
            dncell2_ids[i] = dn2_cell_id

            # Set secondary downstream connection and ratios/demands
            if is_structured:
                split_row, split_col = get_structured_indices(split_cell_idx, nrow, ncol)
                dnID_2d_expand_to_3d[split_row, split_col, 1] = dn2_cell_id

                if is_ibt == 1:
                    # IBT mode: use demand values
                    ibt_demand_3d[split_row, split_col, 0] = 0.0  # Primary gets remainder
                    demand2_val = bif.iloc[i]['demand2']
                    if demand2_val == -9999:
                        demand2_val = 0.0  # Treat missing as 0
                    ibt_demand_3d[split_row, split_col, 1] = demand2_val
                    print(f'  IBT demand: {demand2_val:.2f} m³/s')
                else:
                    # Regular bifurcation: use ratio values
                    ratio1 = bif.iloc[i]['ratio1']
                    ratio2 = bif.iloc[i]['ratio2']

                    if ratio1 == -9999:
                        ratio1 = 0.0
                    if ratio2 == -9999:
                        ratio2 = 0.0

                    # Validate ratios
                    if ratio1 + ratio2 > 1.0001:  # Allow small floating point tolerance
                        raise ValueError(f'Ratios sum to {ratio1 + ratio2:.4f} > 1.0 for row {i}')

                    bifurc_ratio_3d[split_row, split_col, 0] = ratio1
                    bifurc_ratio_3d[split_row, split_col, 1] = ratio2
                    print(f'  Bifurcation ratios: {ratio1:.3f}, {ratio2:.3f}')
            else:
                dnID_1d_expand_to_2d[split_cell_idx, 1] = dn2_cell_id

                if is_ibt == 1:
                    # IBT mode: use demand values
                    ibt_demand_2d[split_cell_idx, 0] = 0.0  # Primary gets remainder
                    demand2_val = bif.iloc[i]['demand2']
                    if demand2_val == -9999:
                        demand2_val = 0.0  # Treat missing as 0
                    ibt_demand_2d[split_cell_idx, 1] = demand2_val
                    print(f'  IBT demand: {demand2_val:.2f} m³/s')
                else:
                    # Regular bifurcation: use ratio values
                    ratio1 = bif.iloc[i]['ratio1']
                    ratio2 = bif.iloc[i]['ratio2']

                    if ratio1 == -9999:
                        ratio1 = 0.0
                    if ratio2 == -9999:
                        ratio2 = 0.0

                    # Validate ratios
                    if ratio1 + ratio2 > 1.0001:  # Allow small floating point tolerance
                        raise ValueError(f'Ratios sum to {ratio1 + ratio2:.4f} > 1.0 for row {i}')

                    bifurc_ratio_2d[split_cell_idx, 0] = ratio1
                    bifurc_ratio_2d[split_cell_idx, 1] = ratio2
                    print(f'  Bifurcation ratios: {ratio1:.3f}, {ratio2:.3f}')

        # Write the updated arrays to NetCDF file
        print('\nWriting updated arrays to NetCDF file...')
        with nc.Dataset(output_file, 'r+') as ds:
            if is_structured:
                ds.variables['bifurc_ratio'][:] = bifurc_ratio_3d
                ds.variables['ibt_demand'][:] = ibt_demand_3d
                ds.variables['dnID'][:] = dnID_2d_expand_to_3d
            else:
                ds.variables['bifurc_ratio'][:] = bifurc_ratio_2d
                ds.variables['ibt_demand'][:] = ibt_demand_2d
                ds.variables['dnID'][:] = dnID_1d_expand_to_2d

        # Fill in the cell IDs in the CSV table
        bif['cell_id'] = cell_ids
        bif['dncell1_id'] = dncell1_ids
        bif['dncell2_id'] = dncell2_ids

        # Write updated CSV file
        output_csv = csv_file.replace('.csv', '_updated.csv')
        bif.to_csv(output_csv, index=False)
        print(f'Updated CSV file written to: {output_csv}')

        # Verify the result
        print('\nVerifying result...')
        verify_bifurc_file(output_file, len(bif), is_structured)


def find_closest_cell(
    lat_array: np.ndarray,
    lon_array: np.ndarray,
    ID_array: np.ndarray,
    target_lat: float,
    target_lon: float,
    is_structured: bool
) -> Tuple[int, int]:
    """
    Find the closest cell to the target latitude and longitude.

    Parameters
    ----------
    lat_array : np.ndarray
        Array of latitudes
    lon_array : np.ndarray
        Array of longitudes
    ID_array : np.ndarray
        Array of cell IDs
    target_lat : float
        Target latitude
    target_lon : float
        Target longitude
    is_structured : bool
        Whether the mesh is structured

    Returns
    -------
    cell_idx : int
        Cell index (-1 if not found)
    cell_id : int
        Cell ID (-1 if not found)
    """
    min_distance = np.inf
    cell_idx = -1
    cell_id = -1

    if is_structured:
        # For structured mesh, search through 2D arrays
        nrow, ncol = lat_array.shape
        for i in range(nrow):
            for j in range(ncol):
                # Calculate distance using simple Euclidean approximation
                lat_diff = lat_array[i, j] - target_lat
                lon_diff = lon_array[i, j] - target_lon
                distance = np.sqrt(lat_diff**2 + lon_diff**2)

                if distance < min_distance:
                    min_distance = distance
                    cell_idx = i * ncol + j  # Convert to 1D index
                    cell_id = int(ID_array[i, j])
    else:
        # For unstructured mesh, search through 1D arrays
        ncells = len(lat_array)
        for i in range(ncells):
            # Calculate distance
            lat_diff = lat_array[i] - target_lat
            lon_diff = lon_array[i] - target_lon
            distance = np.sqrt(lat_diff**2 + lon_diff**2)

            if distance < min_distance:
                min_distance = distance
                cell_idx = i
                cell_id = int(ID_array[i])

    # Check if we found a reasonable match (within ~1 degree)
    if min_distance > 1.0:
        print(f'Warning: Closest cell is {min_distance:.4f} degrees away from target (lat={target_lat:.4f}, lon={target_lon:.4f})')

    print(f'    Found closest cell: distance={min_distance:.4f} degrees, cell_id={cell_id}')

    return cell_idx, cell_id


def get_structured_indices(cell_idx: int, nrow: int, ncol: int) -> Tuple[int, int]:
    """
    Convert 1D cell index to 2D row/col indices for structured mesh.

    Parameters
    ----------
    cell_idx : int
        1D cell index
    nrow : int
        Number of rows
    ncol : int
        Number of columns

    Returns
    -------
    row : int
        Row index
    col : int
        Column index
    """
    col = cell_idx % ncol
    row = cell_idx // ncol
    return row, col


def verify_bifurc_file(filename: str, num_bifurc_points: int, is_structured: bool) -> None:
    """
    Verify the created file has correct structure and data.

    Parameters
    ----------
    filename : str
        Path to NetCDF file to verify
    num_bifurc_points : int
        Expected number of bifurcation points
    is_structured : bool
        Whether the mesh is structured
    """
    print(f'Verification of {filename}:')

    with nc.Dataset(filename, 'r') as ds:
        # Read dimensions
        if 'downstream' in ds.dimensions:
            print(f'  downstream dimension: {ds.dimensions["downstream"].size}')

        # Read dnID and check structure
        dnID_data = ds.variables['dnID'][:]
        bifurc_ratio = ds.variables['bifurc_ratio'][:]
        ibt_demand = ds.variables['ibt_demand'][:]

        if is_structured:
            # Structured mesh verification
            nlon, nlat, ndownstream = dnID_data.shape
            print(f'  Structured mesh - dnID dimensions: {nlon} x {nlat} x {ndownstream}')

            # Count non-missing values in each layer
            for k in range(ndownstream):
                layer = dnID_data[:, :, k]
                valid_count = np.sum((layer > 0) & (layer != -9999))
                print(f'  Layer {k}: {valid_count} valid connections')

            # Count bifurcation points (cells with secondary downstream connections)
            layer2 = dnID_data[:, :, 1]
            bifurc_count = np.sum((layer2 > 0) & (layer2 != -9999))
            print(f'  Bifurcation points found: {bifurc_count}')
        else:
            # Unstructured mesh verification
            ncells, ndownstream = dnID_data.shape
            print(f'  Unstructured mesh - dnID dimensions: {ncells} x {ndownstream}')

            # Count non-missing values in each column
            for k in range(ndownstream):
                column = dnID_data[:, k]
                valid_count = np.sum((column > 0) & (column != -9999))
                print(f'  Column {k}: {valid_count} valid connections')

            # Count bifurcation points (cells with secondary downstream connections)
            column2 = dnID_data[:, 1]
            bifurc_count = np.sum((column2 > 0) & (column2 != -9999))
            print(f'  Bifurcation points found: {bifurc_count}')

        # Check if we have the expected number of bifurcation points
        if bifurc_count == num_bifurc_points:
            print(f'  ✓ Expected {num_bifurc_points} bifurcation points, found {bifurc_count}')
        else:
            print(f'  ✗ Expected {num_bifurc_points} bifurcation points, found {bifurc_count}')

        # Count non-zero ratio and demand values
        nonzero_ratios = np.sum(bifurc_ratio > 0)
        nonzero_demands = np.sum(ibt_demand > 0)
        print(f'  Non-zero bifurc_ratio values: {nonzero_ratios}')
        print(f'  Non-zero ibt_demand values: {nonzero_demands}')

        # Summary of variables created
        print('  Variables created:')
        print('    - dnID: Downstream connectivity with bifurcation support')
        print('    - bifurc_ratio: Split ratios (0.0-1.0, zeros for most cells)')
        print('    - ibt_demand: IBT demands (m³/s, zeros for most cells)')
        print('    - dnID_original: Backup of original dnID')

    print('Verification complete.')


if __name__ == '__main__':
    # Example usage (commented out - uncomment and modify paths as needed)

    # Example 1: Global half-degree mesh with US IBT data
    # input_file = '/compyfs/zhou014/datasets/E3SM_inputs/MOSART_Global_half_20210422.nc'
    # output_file = '/compyfs/zhou014/datasets/E3SM_inputs/MOSART_Global_half_20210422_bifurc_US.nc'
    # csv_file = '/compyfs/zhou014/E3SMv3/bifurcation/MOSART_input_generation/US_IBT_2d.csv'

    # Example 2: SAG mesh with DRB IBT data
    # input_file = '/compyfs/liao313/04model/pyhexwatershed/sag/pyhexwatershed20250708002/hexwatershed/mosart_sag_parameter.nc'
    # csv_file = '/qfs/people/zhou014/E3SMv3/code/MOSART_bifurcation/components/mosart/src/DRB_IBT.csv'
    # output_file = '/compyfs/zhou014/ICoM/dataset/MOSART_sag_parameter_MPAS_bif.nc'

    # create_bifurc_netcdf(input_file, output_file, csv_file)

    pass
