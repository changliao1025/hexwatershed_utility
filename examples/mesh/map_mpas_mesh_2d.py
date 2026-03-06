import os
from osgeo import osr
from hexwatershed_utility.preprocess.mesh.map_mpas_mesh_file import map_mpas_mesh_file
from pyearthviz.map import RasterTileServer

sFilename_mpas_mesh_in = os.path.join('/compyfs/liao313/04model/pyhexwatershed/global/pyflowline20260203001/jigsaw/out','base_mesh.nc')

sFilename_png_out = os.path.join('/qfs/people/liao313/workspace/python/hexwatershed_utility/figures', 'mpas_mesh_ocn6_18coast6lnd12riv6_2d5.png')

#use
dLongitude_center_in = -121
dLatitude_center_in = 31
buffer = 3
aExtent_in = [dLongitude_center_in - buffer * 1.2, dLongitude_center_in + buffer*1.2, dLatitude_center_in - buffer, dLatitude_center_in + buffer]
providers = RasterTileServer.get_available_providers()
print(providers)
map_mpas_mesh_file(sFilename_mpas_mesh_in,
                   sFilename_output_in = sFilename_png_out,
                   iFlag_zebra_in = 1,
              aBasemap_provider_in = ['Esri.NatGeo', 'Esri.Hydro'],
                 iDPI_in = 300,
                 iBasemap_zoom_level_in = None,  # Higher zoom level for better quality (was 11)
                  dLongitude_center_in = dLongitude_center_in,
                    dLatitude_center_in = dLatitude_center_in,
                    aExtent_in = aExtent_in  )