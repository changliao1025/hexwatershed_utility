"""
Inject DEM elevations into MPAS computational mesh NetCDF files.
"""

import math
import os
import netCDF4 as nc
import numpy as np
from osgeo import gdal, ogr, osr


def inject_dem_into_mpas_mesh_netcdf(
    sFilename_mesh_netcdf: str,
    sFilename_mesh_with_dem_netcdf: str,
    aFilename_dem_geotiff: list,
    aFilename_dem_with_ice: list = None,
    sFilename_land_ocean_mask: str = None,
):
    """
    Inject digital elevation model (DEM) values into an MPAS NetCDF mesh.

    Parameters
    ----------
    sFilename_mesh_netcdf : str
        Input MPAS base mesh NetCDF file.
    sFilename_mesh_with_dem_netcdf : str
        Output MPAS mesh NetCDF file with injected elevation fields.
    aFilename_dem_geotiff : list of str
        List of DEM GeoTIFF files to sample from.
    aFilename_dem_with_ice : list of str, optional
        Optional list of DEM with ice GeoTIFF files.
    sFilename_land_ocean_mask : str, optional
        Optional GeoJSON mask file to clip DEMs.
    """
    if not os.path.exists(sFilename_mesh_netcdf):
        raise FileNotFoundError(f"The mesh file does not exist: {sFilename_mesh_netcdf}")

    for sFilename in aFilename_dem_geotiff:
        if not os.path.exists(sFilename):
            raise FileNotFoundError(f"The DEM file does not exist: {sFilename}")

    if aFilename_dem_with_ice:
        for sFilename in aFilename_dem_with_ice:
            if not os.path.exists(sFilename):
                raise FileNotFoundError(f"The DEM with ice file does not exist: {sFilename}")

    if os.path.exists(sFilename_mesh_with_dem_netcdf):
        os.remove(sFilename_mesh_with_dem_netcdf)

    from pyearth.gis.gdal.read.raster.gdal_read_geotiff_file import gdal_read_geotiff_file
    from pyearth.gis.geometry.convert_longitude_range import convert_360_to_180_np

    with nc.Dataset(sFilename_mesh_netcdf, "r") as pDatasets_mesh:
        format = pDatasets_mesh.file_format

        with nc.Dataset(sFilename_mesh_with_dem_netcdf, "w", format=format) as pDatasets_out:
            # Copy dimensions
            for sKey, iValue in pDatasets_mesh.dimensions.items():
                dummy = len(iValue)
                pDatasets_out.createDimension(sKey, dummy)

            # Copy existing variables
            for sKey, aValue in pDatasets_mesh.variables.items():
                outVar = pDatasets_out.createVariable(
                    sKey, aValue.datatype, aValue.dimensions
                )
                outVar[:] = aValue[:]
                outVar.setncatts({k: aValue.getncattr(k) for k in aValue.ncattrs()})

            lonCell0 = pDatasets_mesh.variables["lonCell"]
            latCell0 = pDatasets_mesh.variables["latCell"]

            aLatitudeCell = latCell0[:] / math.pi * 180
            aLongitudeCell = lonCell0[:] / math.pi * 180
            aLongitudeCell_180 = convert_360_to_180_np(aLongitudeCell)
            ncell = len(aLongitudeCell)

            aElevation = np.full(ncell, -9999.0, dtype=float)

            for dem_file in aFilename_dem_geotiff:
                dem_data, pSpatialRef, pGeoTransform = gdal_read_geotiff_file(dem_file)
                # Sample DEM for cells within raster bounds
                x_min = pGeoTransform[0]
                y_max = pGeoTransform[3]
                px_w = pGeoTransform[1]
                px_h = pGeoTransform[5]
                x_max = x_min + px_w * dem_data.shape[1]
                y_min = y_max + px_h * dem_data.shape[0]

                in_bounds = (
                    (aLongitudeCell_180 >= x_min)
                    & (aLongitudeCell_180 <= x_max)
                    & (aLatitudeCell >= y_min)
                    & (aLatitudeCell <= y_max)
                )

                indices = np.where(in_bounds)[0]
                for idx in indices:
                    col = int((aLongitudeCell_180[idx] - x_min) / px_w)
                    row = int((aLatitudeCell[idx] - y_max) / px_h)
                    if 0 <= row < dem_data.shape[0] and 0 <= col < dem_data.shape[1]:
                        val = dem_data[row, col]
                        if val != -9999:
                            aElevation[idx] = val

            elev_var = pDatasets_out.createVariable("elevation", "f8", ("nCells",), fill_value=-9999.0)
            elev_var[:] = aElevation
            elev_var.setncatts({"units": "m", "long_name": "surface elevation"})

    return sFilename_mesh_with_dem_netcdf
