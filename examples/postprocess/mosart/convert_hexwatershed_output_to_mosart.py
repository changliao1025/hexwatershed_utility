


if __name__ == '__main__':
    sRegion = 'sag'

    if sRegion == 'sag':

        sFilename_json_in = '/compyfs/liao313/04model/pyhexwatershed/sag/pyhexwatershed20250708002/hexwatershed/hexwatershed.json'
        # sFilename_mpas_in='/people/liao313/workspace/python/pyhexwatershed_icom/data/sag/input/lnd_mesh.nc'
        sFilename_mosart_parameter_in = '/compyfs/inputdata/rof/mosart/MOSART_Global_half_20210616.nc'
        sFilename_mosart_parameter_out = 'mosart_sag_parameter.nc'
        sFilename_mosart_domain_out = 'mosart_sag_domain.nc'
    else:
        if sRegion == 'susquehanna':
            sFilename_json_in = '/compyfs/liao313/04model/pyhexwatershed/susquehanna/pyhexwatershed20221115001/hexwatershed/hexwatershed.json'
            # sFilename_mpas_in='/qfs/people/liao313/workspace/python/pyhexwatershed_icom/data/susquehanna/input/lnd_cull_mesh.nc'
            sFilename_mosart_parameter_in = '/compyfs/inputdata/rof/mosart/MOSART_Global_half_20210616.nc'
            sFilename_mosart_parameter_out = 'mosart_susquehanna_parameter.nc'
            sFilename_mosart_domain_out = 'mosart_susquehanna_domain.nc'
        else:
            sFilename_json_in = '/compyfs/liao313/04model/pyhexwatershed/columbia/pyhexwatershed20221115003/hexwatershed/hexwatershed.json'
            # sFilename_mpas_in='/compyfs/liao313/00raw/mesh/global/lnd_mesh.nc'
            sFilename_mosart_parameter_in = '/compyfs/inputdata/rof/mosart/MOSART_Global_half_20210616.nc'
            sFilename_mosart_parameter_out = 'mosart_columbia_parameter.nc'
            sFilename_mosart_domain_out = 'mosart_columbia_domain.nc'

    convert_hexwatershed_json_to_mosart_netcdf(sFilename_json_in,
                                               # sFilename_mpas_in, \
                                               sFilename_mosart_parameter_in,
                                               sFilename_mosart_parameter_out,
                                               sFilename_mosart_domain_out)
