import os
import glob
from osgeo import gdal, ogr
gdal.UseExceptions()

from pyearth.system.define_global_variables import *
from pyearth.gis.gdal.read.vector.gdal_get_vector_extent import gdal_get_vector_extent
from pyearth.gis.location.get_hydrosheds_continent_from_extent import get_hydrosheds_continent_from_extent

def find_minimal_hydrobasins_watershed_boundary(sFilename_river_network_in, sFolder_watershed_boundary_in, iFlag_nested_in= False):
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



    # Main function logic
    print(f"Reading river network from: {sFilename_river_network_in}")

    # Get river network bounds and geometry
    aNetwork_bounds, pGeometry_network = gdal_get_vector_extent(sFilename_river_network_in, iFlag_return_union_geometry=1)
    if aNetwork_bounds is None:
        print("Failed to read river network")
        return None

    #get containing folder of the river network file
    sFolder_river_network = os.path.dirname(sFilename_river_network_in)

    print(f"River network bounds: {aNetwork_bounds}")

    # use regex to find the folder what has the name of the region such as 'na' for north america

    sRegion = get_hydrosheds_continent_from_extent(aNetwork_bounds)
    sRegax = '*_' + sRegion + '_*'

    aFolder = glob.glob(os.path.join(sFolder_watershed_boundary_in, sRegax))
    if not aFolder:
        print(f"No watershed folders found for region '{sRegion}' in {sFolder_watershed_boundary_in}")
        return None

    for sFolder in aFolder:
        if not os.path.isdir(sFolder):
            print(f"Skipping non-directory: {sFolder}")
            continue
        else:
            print(f"Searching in watershed folder: {sFolder}")
            break

    iLevel_start = 10
    iLevel_end = 1

    aWkt = list()

    #open river network file
    pDataset_network = ogr.Open(sFilename_river_network_in)
    if pDataset_network is None:
        print(f"Could not open river network file: {sFilename_river_network_in}")
        return None
    pLayer_network = pDataset_network.GetLayer(0)
    if pLayer_network is None:
        print(f"Could not get layer from river network file: {sFilename_river_network_in}")
        return None

    if iFlag_nested_in:
        nSubbasin = pLayer_network.GetFeatureCount() #we assume each feature is a subbasin because they are simplified already by the previous process
        if nSubbasin == 0:
            print(f"No features found in river network file: {sFilename_river_network_in}")
            return None

        for pFeature_network in pLayer_network:

            sWkt = None
            pGeometry_network = pFeature_network.GetGeometryRef()
            if pGeometry_network is None:
                print("Feature has no valid geometry, skipping")
                continue

            #now get the bounding box of the geometry
            aBounds = pGeometry_network.GetEnvelope()
            if aBounds is None:
                print("Failed to get bounding box of network geometry")
                continue

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
                    if pGeometry_watershed.Contains(pGeometry_network):
                        iFound = 1
                        break
                    else:
                        #check whether the river network geometry intersects with the watershed geometry
                        if pGeometry_watershed.Intersects(pGeometry_network):
                            #skip to the next level
                            print(f"River network intersects with watershed at level {sLevel}, skipping to next level")
                            break

                #export the watershed boundary if found
                if iFound == 1:
                    sWkt = pGeometry_watershed.ExportToWkt()
                    aWkt.append(sWkt)
                    #save to a file
                    sFilename_watershed_out = os.path.join(sFolder_river_network, f"subbasin_{pFeature_network.GetFID()}_watershed_boundary_{sLevel}.geojson")
                    pDriver = ogr.GetDriverByName('GeoJSON')
                    if os.path.exists(sFilename_watershed_out):
                        os.remove(sFilename_watershed_out)
                    pDataset_out = pDriver.CreateDataSource(sFilename_watershed_out)
                    if pDataset_out is None:
                        print(f"Could not create output dataset: {sFilename_watershed_out}")
                        continue
                    pLayer_out = pDataset_out.CreateLayer("watershed_boundary", geom_type=ogr
                    .wkbPolygon)
                    if pLayer_out is None:
                        print(f"Could not create layer in output dataset: {sFilename_watershed_out}")
                        continue
                    pFeature_out = ogr.Feature(pLayer_out.GetLayerDefn())
                    pGeometry_out = ogr.CreateGeometryFromWkt(sWkt)
                    if pGeometry_out is None:
                        print(f"Could not create geometry from WKT: {sWkt}")
                        continue
                    pFeature_out.SetGeometry(pGeometry_out)
                    if pLayer_out.CreateFeature(pFeature_out) != ogr.OGRERR_NONE:
                        print(f"Could not create feature in output layer: {sFilename_watershed_out}")
                        continue
                    pFeature_out = None
                    pDataset_out = None
                    print(f"Watershed boundary for subbasin {pFeature_network.GetFID()} at level {sLevel} saved to {sFilename_watershed_out}")
                    break
                else:
                    continue



        #merge all the wkt into one geometry
        if len(aWkt) == 0:
            print("No watershed boundaries found for any subbasin")
            return None
        else:
            pGeometry_merged = None
            for sWkt in aWkt:
                pGeometry = ogr.CreateGeometryFromWkt(sWkt)
                if pGeometry_merged is None:
                    pGeometry_merged = pGeometry
                else:
                    pGeometry_merged = pGeometry_merged.Union(pGeometry)
            if pGeometry_merged is None:
                print("Failed to merge watershed geometries")
                return None
            sWkt_merged = pGeometry_merged.ExportToWkt()
            return sWkt_merged
    else:
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
                if pGeometry_watershed.Contains(pGeometry_network):
                    iFound = 1
                    break
                else:
                    #check whether the river network geometry intersects with the watershed geometry
                    if pGeometry_watershed.Intersects(pGeometry_network):
                        #skip to the next level
                        print(f"River network intersects with watershed at level {sLevel}, skipping to next level")
                        break
            #export the watershed boundary if found
            if iFound == 1:
                sWkt = pGeometry_watershed.ExportToWkt()
                aWkt.append(sWkt)
                break
            if sWkt is None:
                print("No containing watershed boundary found for", pGeometry_network.ExportToWkt())
                return None
            else:
                #print(f"Found containing watershed boundary at level {sLevel}: {sWkt}")
                return sWkt
