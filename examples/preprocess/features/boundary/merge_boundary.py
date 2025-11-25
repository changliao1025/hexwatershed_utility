import os, sys
import numpy as np
from osgeo import gdal, ogr, osr

sPath_project = '/qfs/people/liao313/workspace/python/hexwatershed_utility/'
#add the project path of the pythonpath
sys.path.append(sPath_project)
from hexwatershed_utility.codes.merge_features import merge_features

sFilename_polygon_in='/qfs/people/liao313/workspace/python/hexwatershed_utility/data/conus/US_State_Boundaries_CONUS.geojson'
sExtension_vector = '.geojson'
#merge polygon as one

print('The polygon contains more than one polygon, the program will attempt to merge them as one!')
            #obtain the file extension
sFilename_clip_new = '/qfs/people/liao313/workspace/python/hexwatershed_utility/data/conus/conus_boundary_new.geojson'
merge_features(sFilename_polygon_in, sFilename_clip_new)