import os
import glob

from osgeo import gdal, ogr, osr
gdal.UseExceptions()

from pyearth.system.define_global_variables import *

def find_minimal_watershed_boundary_hydrobasin(sRegion, sFilename_river_network_in, sFolder_watershed_boundary_in):
    """
    Find the minimal watershed boundary that completely contains a river network.

    Searches from low level (small watersheds) to high level (large watersheds)
    to find the smallest boundary that encompasses the entire river network.

    Args:
        sFilename_river_network_in (str): Path to river network GeoJSON file
        sFolder_watershed_boundary_in (str): Path to folder containing watershed boundary shapefiles
                                           organized by levels (e.g., level_01, level_02, etc.)

    Returns:
        dict: Dictionary containing:
            - 'level': The watershed level (e.g., 'level_03')
            - 'filename': Path to the shapefile containing the boundary
            - 'feature_id': ID/index of the feature that contains the network
            - 'geometry': The boundary geometry as OGR geometry object
            - 'bounds': Bounding box of the containing watershed (minx, miny, maxx, maxy)
            - 'area': Area of the containing watershed
    """

    def get_river_network_bounds(river_file):
        """Extract the bounding box and union geometry of the entire river network"""
        try:
            # Open river network file
            dataset = ogr.Open(river_file)
            if dataset is None:
                print(f"Could not open {river_file}")
                return None, None

            layer = dataset.GetLayer(0)
            if layer is None:
                print("Could not get layer from river network file")
                return None, None

            # Create a geometry collection to hold all river segments
            network_union = ogr.Geometry(ogr.wkbMultiLineString)

            # Get extent of all features
            minx, maxx, miny, maxy = layer.GetExtent()
            network_bounds = (minx, miny, maxx, maxy)

            # Collect all geometries
            layer.ResetReading()
            for feature in layer:
                geometry = feature.GetGeometryRef()
                if geometry is not None:
                    # Clone geometry to avoid memory issues
                    geom_clone = geometry.Clone()

                    # Handle different geometry types
                    geom_type = geom_clone.GetGeometryType()

                    if geom_type == ogr.wkbLineString:
                        network_union.AddGeometry(geom_clone)
                    elif geom_type == ogr.wkbMultiLineString:
                        for i in range(geom_clone.GetGeometryCount()):
                            line = geom_clone.GetGeometryRef(i)
                            network_union.AddGeometry(line.Clone())

            # Close dataset
            dataset = None

            return network_bounds, network_union

        except Exception as e:
            print(f"Error reading river network: {e}")
            return None, None







    # Main function logic
    print(f"Reading river network from: {sFilename_river_network_in}")

    # Get river network bounds and geometry
    network_bounds, network_geometry = get_river_network_bounds(sFilename_river_network_in)
    if network_bounds is None:
        print("Failed to read river network")
        return None

    print(f"River network bounds: {network_bounds}")

    # use regax to find the folder what has the name of the region such as 'na' for north america

    sRegax = '*' + sRegion + '*'

    aFolder = glob.glob(os.path.join(sFolder_watershed_boundary_in, sRegax))
    if not aFolder:
        print(f"No watershed folders found for region '{sRegion}' in {sFolder_watershed_boundary_in}")
        return None

    for sFolder in aFolder:
        if not os.path.isdir(sFolder):
            print(f"Skipping non-directory: {sFolder}")
            continue

    iLevel_start = 4
    iLevel_end = 1
    sWkt = None
    for iLevel in range(iLevel_start, iLevel_end, -1):
        sLevel = f"lev{iLevel:02d}"
        sFilename_reg = 'hybas_lake_*' + sLevel + '*_*.shp'
        aFilename = glob.glob(os.path.join(sFolder, sFilename_reg))
        if not aFilename:
            print(f"No shapefiles found for level {sLevel} in {sFolder}")
            continue
        sFilename_full = aFilename[0]  # Use the first matching shapefile
        #sFilename_full = os.path.join(sFolder, sFilename)

        iFound = 0
        #use the file to search the geometry
        pDataset_watershed = ogr.Open(sFilename_full)
        if pDataset_watershed is None:
            print(f"Could not open watershed file: {sFilename_full}")
            continue
        pLayer_watershed = pDataset_watershed.GetLayer(0)
        if pLayer_watershed is None:
            print(f"Could not get layer from watershed file: {sFilename_full}")
            continue
        # Iterate through watershed features
        pLayer_watershed.ResetReading()
        for pFeature_watershed in pLayer_watershed:
            pGeometry_watershed = pFeature_watershed.GetGeometryRef()
            if pGeometry_watershed is None:
                continue

            # Check if the watershed geometry contains the river network geometry
            if pGeometry_watershed.Contains(network_geometry):
                iFound = 1
                break
            else:
                #check whether the river network geometry intersects with the watershed geometry
                if pGeometry_watershed.Intersects(network_geometry):
                    #skip to the next level
                    print(f"River network intersects with watershed at level {sLevel}, skipping to next level")
                    break

        #export the watershed boundary if found
        if iFound == 1:
            sWkt = pGeometry_watershed.ExportToWkt()
            break



    return sWkt



if __name__ == "__main__":
    sRegion = 'na'
    sFolder_watershed_boundary_in = '/compyfs/liao313/00raw/hydrology/hydrosheds/hydrobasin'
    sFolder_out  = '/compyfs/liao313/04model/pyhexwatershed/northamerica/watershed_boundary'
    sFolder_in = '/compyfs/liao313/04model/pyhexwatershed/northamerica/river_network'
    pDriver = ogr.GetDriverByName('GeoJSON')

    aFilename = list()
    aFilename.append(os.path.join(sFolder_in, 'delaware.geojson'))
    aFilename.append(os.path.join(sFolder_in, 'susquehanna.geojson'))


    #for i in range(1, 11):
    for sFilename_river_network_in in aFilename:
        #sBasin  = f'{i:04d}'
        #sFilename_river_network_in = '/compyfs/liao313/04model/pyhexwatershed/northamerica/river_network/HydroRIVERS_v10_na_simplified_1.25E+04_3.12E+09_'+sBasin+'_outlet_simplified.geojson'

        wkt = find_minimal_watershed_boundary_hydrobasin(sRegion, sFilename_river_network_in, sFolder_watershed_boundary_in)

        #save a geojson file
        #sFilename_out = os.path.join(sFolder_out, f"watershed_boundary_{sBasin}.geojson")
        sBasin = os.path.splitext(os.path.basename(sFilename_river_network_in))[0]
        sFilename_out = os.path.join(sFolder_out, f"watershed_boundary_{sBasin}.geojson")
        if os.path.exists(sFilename_out):
            os.remove(sFilename_out)
        pDataset_out = pDriver.CreateDataSource(sFilename_out)
        if pDataset_out is None:
            print(f"Could not create output dataset: {sFilename_out}")
            continue
        pLayer_out = pDataset_out.CreateLayer("watershed_boundary", geom_type=ogr.wkbPolygon)
        if pLayer_out is None:
            print(f"Could not create layer in output dataset: {sFilename_out}")
            continue
        pFeature_out = ogr.Feature(pLayer_out.GetLayerDefn())
        pGeometry_out = ogr.CreateGeometryFromWkt(wkt)
        if pGeometry_out is None:
            print(f"Could not create geometry from WKT: {wkt}")
            continue
        pFeature_out.SetGeometry(pGeometry_out)
        if pLayer_out.CreateFeature(pFeature_out) != ogr.OGRERR_NONE:
            print(f"Could not create feature in output layer: {sFilename_out}")
            continue
        pFeature_out = None
        pDataset_out = None

        print(f"Minimal watershed boundary for basin {sBasin} saved to {sFilename_out}")

    print(f"Minimal watershed boundary WKT: {wkt}")

