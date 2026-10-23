import json
import os
from pathlib import Path

import netCDF4 as nc
import numpy as np
from pyearth.gis.geometry.calculate_polygon_area import calculate_polygon_area
from pyesm.mesh.unstructured.esm_convert_unstructured_domain_file_to_scripgrid_file import (
    esm_convert_unstructured_domain_file_to_scripgrid_file,
)
from pyesm.mesh.unstructured.esm_create_unstructured_domain_file_full import (
    esm_create_unstructured_domain_file_full,
)


def convert_hexwatershed_json_to_land_domain_file(
    sFilename_json_in, sFilename_land_domain_out
):
    aLongitude = list()
    aLatitude = list()
    aLonV_region = list()
    aLatV_region = list()
    aCellID = list()
    aArea = list()
    with open(sFilename_json_in) as json_file:
        data = json.load(json_file)
        ncell = len(data)
        lID = 1
        for i in range(ncell):
            pcell = data[i]
            lID = lID + 1
            # aArea.append(dArea)
            aLongitude.append(float(pcell["dLongitude_center_degree"]))
            aLatitude.append(float(pcell["dLatitude_center_degree"]))
            aCellID.append(int(pcell["lCellID"]))

            dummy_vertex = pcell["vVertex"]
            aVertex_lon = np.full(9, -9999, float)
            aVertex_lat = np.full(9, -9999, float)
            nVertex = len(dummy_vertex)
            for j in range(nVertex):
                aVertex_lon[j] = dummy_vertex[j]["dLongitude_degree"]
                aVertex_lat[j] = dummy_vertex[j]["dLatitude_degree"]

            aLonV_region.append(aVertex_lon)
            aLatV_region.append(aVertex_lat)
            dArea = float(pcell["dArea"])
            aArea.append(dArea)

        pass

    aLongitude = np.array(aLongitude).reshape(ncell)
    aLatitude = np.array(aLatitude).reshape(ncell)
    nCell = len(aCellID)
    aCellID = np.array(aCellID).reshape(ncell)
    aArea = np.array(aArea).reshape(ncell)

    aLonV_region = np.array(aLonV_region)  # .T
    aLatV_region = np.array(aLatV_region)  # .T

    nrow = nCell
    aLongitude.shape = (nrow, 1)
    aLatitude.shape = (nrow, 1)

    nrow, nvertex = aLonV_region.shape
    aLonV_region.shape = (nrow, 1, nvertex)
    aLatV_region.shape = (nrow, 1, nvertex)

    esm_create_unstructured_domain_file_full(
        aLongitude,
        aLatitude,
        aLonV_region,
        aLatV_region,
        sFilename_land_domain_out,
        aArea_in=aArea,
    )

    esm_convert_unstructured_domain_file_to_scripgrid_file(
        sFilename_land_domain_out,
        sFilename_land_domain_out.replace(".nc", ".scripgrid.nc"),
    )

    return
