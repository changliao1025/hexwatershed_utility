import os
from pathlib import Path
import json
import numpy as np
import netCDF4 as nc
from pye3sm.mesh.unstructured.e3sm_create_unstructured_domain_file_full import e3sm_create_unstructured_domain_file_full
from pye3sm.mesh.unstructured.e3sm_convert_unstructured_domain_file_to_scripgrid_file import e3sm_convert_unstructured_domain_file_to_scripgrid_file

#because elm domain only need to cell center and vertex info, it is best to use the mesh info file

def convert_pyflowline_mesh_to_elm_domain_file(sFilename_mesh_info_in, sFilename_domain_out):
    netcdf_format = 'NETCDF4'

    if os.path.exists(sFilename_domain_out):
        os.remove(sFilename_domain_out)

    pDatasets_out = nc.Dataset(
        sFilename_domain_out, "w", format=netcdf_format)

    aID = list()
    aLonV_region=list()
    aLatV_region=list()
    aLongitude    =list()
    aLatitude    =list()
    aArea   =list()
    aCellID = list()
    lID = 1
    with open(sFilename_mesh_info_in) as json_file:
        data = json.load(json_file)
        ncell = len(data)
        nrow = ncell
        for i in range(ncell):
            pcell = data[i]
            aID.append(lID)
            lID = lID + 1

            #aArea.append(dArea)
            aLongitude.append(float(pcell['dLongitude_center_degree']))
            aLatitude.append(float(pcell['dLatitude_center_degree']))
            aCellID.append(int(pcell['lCellID']))

            dummy_vertex = pcell['aVertex']
            aVertex_lon = np.full(9, -9999, float)
            aVertex_lat = np.full(9, -9999, float)
            nVertex = len(dummy_vertex)
            for j in range(nVertex):
                aVertex_lon[j] = dummy_vertex[j]['dLongitude_degree']
                aVertex_lat[j] = dummy_vertex[j]['dLatitude_degree']

            aLonV_region.append(aVertex_lon)
            aLatV_region.append(aVertex_lat)
            dArea = float(pcell['dArea'])
            aArea.append(dArea)

    #convert to numpy array
    aLongitude = np.array(aLongitude)
    aLatitude = np.array(aLatitude)
    aLonV_region = np.array(aLonV_region)
    aLatV_region = np.array(aLatV_region)
    aArea = np.array(aArea)
    aLongitude.shape = (nrow, 1)
    aLatitude.shape = (nrow, 1)

    nrow, nvertex = aLonV_region.shape
    aLonV_region.shape = (nrow, 1, nvertex)
    aLatV_region.shape = (nrow, 1, nvertex)

    e3sm_create_unstructured_domain_file_full(aLongitude, aLatitude,  aLonV_region, aLatV_region,
                                              sFilename_domain_out, aArea_in=aArea)

    e3sm_convert_unstructured_domain_file_to_scripgrid_file(sFilename_domain_out, sFilename_domain_out.replace('.nc', '.scripgrid.nc'))

