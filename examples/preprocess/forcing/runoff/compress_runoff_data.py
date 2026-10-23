import numpy as np
from netCDF4 import Dataset  # read netcdf
from osgeo import gdal, osr  # the default operator
from pyearth.gis.gdal.write.gdal_write_envi_file import (
    gdal_write_envi_file_multiple_band,
)
from pyearth.gis.gdal.write.gdal_write_geotiff_file import (
    gdal_write_geotiff_file_multiple_band,
)
from pyearth.system.define_global_variables import *
from scipy.io import loadmat

sWorkspace_runoff = "/qfs/people/xudo627/ming/runoff"
sWorkspace_runoff_out = "/compyfs/liao313/00raw/"


def compress_runoff_data():
    ncolumn = 720
    nrow = 300
    dLon_min = -179.75
    dLat_max = 89.75
    dResolution_x = 0.5
    dResolution_y = 0.5
    pHeaderParameters = {}
    pHeaderParameters["ncolumn"] = f"{ncolumn:0d}"
    pHeaderParameters["nrow"] = f"{nrow:0d}"
    pHeaderParameters["ULlon"] = f"{dLon_min - 0.5 * dResolution_x:0f}"
    pHeaderParameters["ULlat"] = f"{dLat_max + 0.5 * dResolution_y:0f}"
    pHeaderParameters["pixelSize"] = f"{dResolution_x:0f}"
    pHeaderParameters["nband"] = "1"
    pHeaderParameters["offset"] = "0"
    pHeaderParameters["data_type"] = "4"
    pHeaderParameters["bsq"] = "bsq"
    pHeaderParameters["byte_order"] = "0"
    pHeaderParameters["missing_value"] = "-9999"
    iYear_end = 2009
    iYear_start = 1979
    nday = (iYear_end - iYear_start + 1) * 365
    aGrid_stack = np.full((nday, nrow, ncolumn), -9999.0, dtype=float)
    k = 0
    for i in range(1979, 2009 + 1, 1):
        si = f"{i:04d}"
        for j in range(1, 366, 1):
            sj = f"{j:0d}"

            sFilename = sWorkspace_runoff + "/" + "RUNOFF05_" + si + "_" + sj + ".mat"
            print(sFilename)
            dummy_data = loadmat(sFilename)
            dummy = np.transpose(dummy_data["ro05"])
            nan_index = np.where(np.isfinite(dummy) == False)
            dummy[nan_index] = -9999
            aGrid_stack[k, :, :] = dummy
            k = k + 1

    pSpatial = osr.SpatialReference()
    pSpatial.ImportFromEPSG(4326)
    sFilename_envi = sWorkspace_runoff_out + slash + "runoff" + sExtension_envi

    gdal_write_envi_file_multiple_band(
        sFilename_envi,
        aGrid_stack,
        float(pHeaderParameters["pixelSize"]),
        float(pHeaderParameters["ULlon"]),
        float(pHeaderParameters["ULlat"]),
        -9999.0,
        pSpatial,
    )


if __name__ == "__main__":
    compress_runoff_data()
