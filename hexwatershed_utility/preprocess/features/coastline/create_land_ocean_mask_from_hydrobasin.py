import os
import glob
from osgeo import gdal, ogr
gdal.UseExceptions()

from pyearth.toolbox.management.vector.merge_files import merge_files
from pyearth.toolbox.management.vector.remove_small_polygon import remove_small_polygon
from pyearth.toolbox.management.vector.merge_features import merge_features
from pyearth.toolbox.geometry.create_gcs_buffer_zone import create_buffer_zone_polygon_file
from pyearth.toolbox.conversion.convert_vector_format import convert_vector_format
from pyearth.toolbox.conversion.convert_vector_to_global_raster import convert_vector_to_global_raster
def create_land_ocean_mask_from_hydrobasin(sWorkspace_coastline_output,
                                                                             sWorkspace_watershed_boundary_in,
                                                                             dResolution_x_in, dResolution_y_in,
                                                                             dThreshold_area_island,
                                                                             dResolution_coastline_buffer):

    pDriver_shapefile = ogr.GetDriverByName('ESRI Shapefile')
    if pDriver_shapefile is None:
        print("ESRI Shapefile driver not available.")
        return


    aRegion = list()
    aRegion.append('af')
    aRegion.append('ar')
    aRegion.append('as')
    aRegion.append('au')
    aRegion.append('eu')
    aRegion.append('gr')
    aRegion.append('na')
    aRegion.append('sa')
    aRegion.append('si')

    #example hybas_lake_af_lev01-12_v1c/
    aWatershed_boundary = list()
    for sRegion in aRegion:
        sFolder = 'hybas_lake_' + sRegion + '_lev01-12_v1c'
        sFolder_in = os.path.join(sWorkspace_watershed_boundary_in, sFolder)

        #only used level 02
        sLevel = 'lev03'
        sFilename_reg = 'hybas_lake_*' + sLevel + '*_*.shp'
        aFilename = glob.glob(os.path.join(sFolder_in, sFilename_reg))
        if not aFilename:
            print(f"No shapefiles found for level {sLevel} in {sFolder_in}")
            continue
        sFilename_full = aFilename[0]  # Use the first matching shapefile
        aWatershed_boundary.append(sFilename_full)

    # Merge all the geojson files into a single file
    sFilename_merge = os.path.join(sWorkspace_coastline_output, 'hybas_lake_all_lev03.parquet')
    merge_files(aWatershed_boundary, sFilename_merge, sFormat='Parquet')

    sFilename_wo_island = os.path.join(sWorkspace_coastline_output, 'land_ocean_mask_wo_island.parquet')
    remove_small_polygon(sFilename_merge, sFilename_wo_island, dThreshold_area_island )

    exit()
    sFilename_vector_merged_raw = os.path.join(sWorkspace_coastline_output, 'land_ocean_mask_wo_island_merged_raw.parquet')
    merge_features(sFilename_wo_island, sFilename_vector_merged_raw)

    sFilename_tif_wo_island = os.path.join(sWorkspace_coastline_output, 'land_ocean_mask_wo_island.tif')
    convert_vector_to_global_raster(sFilename_vector_merged_raw,
                                    sFilename_tif_wo_island,
                                    dResolution_x_in,
                                    dResolution_y_in,
                                    iFlag_boundary_only_in = 0,
                                    dFill_value_in = 2)

    exit
    sFilename_parquet_buffer = os.path.join(sWorkspace_coastline_output, 'land_ocean_mask_buffer.parquet')
    create_buffer_zone_polygon_file(sFilename_vector_merged_raw, sFilename_parquet_buffer,
                                          dBuffer_distance_in = dResolution_coastline_buffer * 1000)

    sFilename_vector_parquet_merged = os.path.join(sWorkspace_coastline_output, 'land_ocean_mask_buffer_merged.parquet')
    merge_features(sFilename_parquet_buffer, sFilename_vector_parquet_merged)


    return sFilename_tif_wo_island, sFilename_vector_parquet_merged
