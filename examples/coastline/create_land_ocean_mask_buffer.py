import os, sys
import numpy as np
from osgeo import gdal, ogr, osr
from scipy import ndimage
from pyearth.system.define_global_variables import *

from pyearth.toolbox.data.ocean.define_land_ocean_mask import create_land_ocean_vector_mask
from pyearth.toolbox.conversion.convert_vector_to_geojson import convert_vector_to_geojson
from pyearth.toolbox.conversion.convert_vector_to_global_raster import convert_vector_to_global_raster
from pyearth.toolbox.management.vector.remove_small_polygon import remove_small_polygon
from pyearth.gis.gdal.read.raster.gdal_read_geotiff_file import gdal_read_geotiff_file
from pyearth.toolbox.conversion.vectorize_raster import vectorize_raster
from pyearth.gis.gdal.write.raster.gdal_write_geotiff_file import gdal_write_geotiff_file
from pyearth.toolbox.management.vector.fields import get_field_and_value, add_field_to_vector_file
from pyearth.toolbox.management.vector.merge_features import merge_features
from pyearth.toolbox.geometry.create_gcs_buffer_zone import create_gcs_buffer_zone_polygon
def fill_holes(mask):
    """
    Fills holes in a binary NumPy mask and then finds its edge.

    Args:
        mask (np.ndarray): A 2D binary NumPy array (dtype=bool or 0/1).

    Returns:
        np.ndarray: A boolean NumPy array of the same shape as the input mask,
                    where True indicates an edge pixel and False otherwise.
    """
    if mask.ndim != 2:
        raise ValueError("Input mask must be a 2D array.")

    # Convert to boolean if it's not already
    mask = mask.astype(bool)

    # Fill holes
    filled_mask = ndimage.binary_fill_holes(mask)

    return filled_mask

def find_edge(filled_mask):

    if filled_mask.ndim != 2:
        raise ValueError("Input mask must be a 2D array.")

    # Convert to boolean if it's not already
    filled_mask = filled_mask.astype(bool)

    # Pad the filled mask to handle edges correctly
    padded_mask = np.pad(filled_mask, pad_width=1, mode='constant', constant_values=False)

    # Create a result array initialized to False
    edge = np.zeros_like(filled_mask, dtype=bool)

    # Iterate through the original filled mask (excluding padding)
    for r in range(filled_mask.shape[0]):
        for c in range(filled_mask.shape[1]):
            # Check if the current pixel is part of the object (in the filled mask)
            if padded_mask[r + 1, c + 1]:
                # Check its neighbors in the padded mask
                neighbors = [
                    padded_mask[r, c + 1],     # Top
                    padded_mask[r + 2, c + 1],   # Bottom
                    padded_mask[r + 1, c],     # Left
                    padded_mask[r + 1, c + 2],   # Right
                    # Optional: Include diagonals if you define edges that way
                    padded_mask[r, c],         # Top-Left
                    padded_mask[r, c + 2],       # Top-Right
                    padded_mask[r + 2, c],       # Bottom-Left
                    padded_mask[r + 2, c + 2]   # Bottom-Right
                ]
                # If any neighbor is False (background), then the current pixel is an edge
                if not all(neighbors):
                    edge[r, c] = True

    return edge

sFilename_geojson_raw = '/qfs/people/liao313/data/hexwatershed/global/vector/region.geojson'

sRegion = 'na'  # Define the region, e.g., 'na' for North America
sFilename_geojson_na_raw = '/compyfs/liao313/00raw/hydrology/hydroshed/hydrobasin/hybas_lake_na_lev01-12_v1c/hybas_lake_na_lev01_v1c.shp'

sWorkspace_vector = '/qfs/people/liao313/data/hexwatershed/global/vector/'
sWorkspace_raster = '/qfs/people/liao313/data/hexwatershed/global/raster/'

sWorkspace_output0 = '/qfs/people/liao313/workspace/python/hexwatershed_utility/hexwatershed_utility/coastline/output/'
sWorkspace_output1 = os.path.join(sWorkspace_output0, sRegion)
if not os.path.exists(sWorkspace_output1):
    os.makedirs(sWorkspace_output1)

#step 0, make a copy of the raw vector file in geojson format
sFilename_geojson_na = os.path.join(sWorkspace_output1, sRegion + '_hydrobasin.geojson')
convert_vector_to_geojson(sFilename_geojson_na_raw, sFilename_geojson_na)

#define target resolution, this resolution will used as the JIGSAW spacing parameter
#for example if the resolution is 5km, then

iResolution_km = 4
sResolution = "{:d}".format(iResolution_km) + 'km'

sWorkspace_output = os.path.join(sWorkspace_output1, sResolution)
if not os.path.exists(sWorkspace_output):
    os.makedirs(sWorkspace_output)

dResolution_x_in = 30.0/3600 * iResolution_km
dResolution_y_in = dResolution_x_in

nrow = int(180 / dResolution_y_in)
ncolumn = int(360 / dResolution_x_in)
print(ncolumn, nrow )

#step 1: record attribute from the MPAS tools
aField, aValue = get_field_and_value(sFilename_geojson_raw)

#step 2: remove islands in hydrobasin
iFlag_remove_small_island = 1
dThreshold_area = 3.0E6  #in square km for small islans
sThreshold_island = "{:.2E}".format(dThreshold_area)+'km2'
sTemp = 'land_ocean_mask_wo_island.geojson'
sFilename_geojson_out = os.path.join(sWorkspace_output, sTemp)
remove_small_polygon(sFilename_geojson_na, sFilename_geojson_out, dThreshold_area )

#step 3, add a buffer to polygon using the resolution

sFilename_geojson_buffer = os.path.join(sWorkspace_output, 'land_ocean_mask_wo_island_buffer_' + sResolution + '.geojson')
create_gcs_buffer_zone_polygon(sFilename_geojson_out, sFilename_geojson_buffer,
                                      dBuffer_distance_in = iResolution_km* 1000)


#step 3, convert the vector to raster and remove the holes
sTemp = 'land_ocean_mask_wo_island_' + sResolution +'.tif'
sFilename_tif_out = os.path.join(sWorkspace_output, sTemp)
print('converting to raster', sFilename_tif_out)
convert_vector_to_global_raster(sFilename_geojson_buffer,
                                sFilename_tif_out,
                                 dResolution_x_in, dResolution_y_in,
                                iFlag_boundary_only_in = 1,
                                 dFill_value_in = 1)
#step 4: convert the filled mask back to a new raster file
#read the raster as a 2D numpy array
dum = gdal_read_geotiff_file(sFilename_tif_out)
aData = dum['dataOut']
dPixelWidth_in = dum['pixelWidth']
dPixelHeight_in = dum['pixelHeight']
dOriginX_in = dum['originX']
dOriginY_in = dum['originY']
dMissing_value_in = dum['missingValue']
pProjection_in = dum['projection']
#remove the holes in the 2D numpy array
filled_mask = fill_holes(aData)
sFilename_nohole = os.path.join(sWorkspace_output, 'land_ocean_mask_wo_island_wo_hole_' + sResolution + '.tif')
gdal_write_geotiff_file(sFilename_nohole,
                            filled_mask,
                            dPixelWidth_in,
                            dPixelHeight_in,
                            dOriginX_in,
                            dOriginY_in,
                            dMissing_value_in,
                            pProjection_in,
                            datatype=gdal.GDT_Byte)

#step 5: convert it back to vector
sFilename_vector_out = os.path.join(sWorkspace_output, 'land_ocean_mask_wo_island_wo_hole_'  + sResolution + '.shp')
vectorize_raster(sFilename_nohole, sFilename_vector_out, sFieldname='land', sFieldtype=ogr.OFTInteger)

#step 6: convert it to geojson and merge the features as one feature
sFilename_vector_geojson = os.path.join(sWorkspace_output, 'land_ocean_mask_wo_island_wo_hole_' + sResolution + '.geojson')
convert_vector_to_geojson(sFilename_vector_out, sFilename_vector_geojson)
sFilename_vector_geojson_merged = os.path.join(sWorkspace_output, 'land_ocean_mask_wo_island_wo_hole_merged_' + sResolution + '.geojson')
merge_features(sFilename_vector_geojson, sFilename_vector_geojson_merged)

#step 7: add the field and value back into the geojson, so it can used as the MPAS tools input
add_field_to_vector_file(sFilename_vector_geojson_merged, aField, aValue)

#the merged geojson file is the final land ocean mask that should be used in the MPAS tools
#for jigsaw coastline vector and raster, we dont need to use the merged geojson file
#step 8, generate the jigsaw spacing raster
edge_mask = find_edge(filled_mask)
#create a new raster with the edge mask
aData_new = np.zeros_like(filled_mask, dtype=np.uint8)
aData_new[filled_mask] = 2  # Set land mask to 2
aData_new[edge_mask] = 1  # Set edge mask to 1
sFilename_jigsaw_raster = os.path.join(sWorkspace_output, 'jigsaw_coastline_mask_' + sResolution + '.tif')
gdal_write_geotiff_file(sFilename_jigsaw_raster,
                            aData_new,
                            dPixelWidth_in,
                            dPixelHeight_in,
                            dOriginX_in,
                            dOriginY_in,
                            dMissing_value_in,
                            pProjection_in,
                            datatype=gdal.GDT_Byte)



print('Finished creating land ocean mask in raster.')