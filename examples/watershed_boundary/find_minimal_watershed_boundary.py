import os
import glob
import re
from osgeo import gdal, ogr
gdal.UseExceptions()
from hexwatershed_utility.preprocess.features.watershed_boundary.find_minimal_hydrobasins_watershed_boundary import find_minimal_hydrobasins_watershed_boundary


dResolution_land = 10
dDistance_tolerance_in = dResolution_land * 1.0E3


dDrainage_area_threshold_in = dResolution_land * dResolution_land *10 * 1.0E6 #km2
sDistance_tolerance = "{:.2E}".format(dDistance_tolerance_in)
sDrainage_area_threshold = "{:.2E}".format(dDrainage_area_threshold_in)



sFolder_watershed_boundary_in = '/compyfs/liao313/00raw/hydrology/hydrosheds/hydrobasin'
sFolder_out  = '/compyfs/liao313/04model/pyhexwatershed/global/watershed_boundary'
sFolder_in = '/compyfs/liao313/04model/pyhexwatershed/global/river_network'
pDriver = ogr.GetDriverByName('GeoJSON')
aFilename = list()
#aFilename.append(os.path.join(sFolder_in, 'delaware.geojson'))
#aFilename.append(os.path.join(sFolder_in, 'susquehanna.geojson'))

#find all the geojson files in the folder
aFilename = glob.glob(os.path.join(sFolder_in, '*.geojson'))

for i in range(12, 13):
#for sFilename_river_network_in in aFilename:
    sBasin  = f'{i:04d}'
    #sFilename_river_network_in = '/compyfs/liao313/04model/pyhexwatershed/northamerica/river_network/HydroRIVERS_v10_na_simplified_1.25E+04_3.12E+09_'+sBasin+'_outlet_simplified.geojson'
    #use glob to find the file
    aFile = glob.glob(os.path.join(sFolder_in, f'*{sBasin}.geojson'))
    if len(aFile) == 0:
        print(f"No file found for basin {sBasin}")
        continue
    sFilename_river_network_in = aFile[0]

    wkt = find_minimal_hydrobasins_watershed_boundary(sFilename_river_network_in, sFolder_watershed_boundary_in)
    #save a geojson file
    #sFilename_out = os.path.join(sFolder_out, f"watershed_boundary_{sBasin}.geojson")
    #sBasin = os.path.splitext(os.path.basename(sFilename_river_network_in))[1]
    if wkt is not None:
        #sFilename_out = os.path.join(sFolder_out, f"watershed_boundary_{sBasin}.geojson")
        #use threshold to define the output filename
        sFilename_out = os.path.join(sFolder_out, f"watershed_boundary_{sDistance_tolerance}_{sDrainage_area_threshold}_{sBasin}.geojson")

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

    else:
        sBasin = None

print(f"Minimal watershed boundary WKT: {wkt}")