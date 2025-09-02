import numpy as np
from pyearth.gis.geometry.split_polygon_cross_idl import split_polygon_cross_idl
aCoords_gcs = [[-179.93378563 ,  71.07872049], [-179.91974662 ,  71.14487311], [ 179.96781872 ,  71.17125621], [ 179.78364363 ,  71.1494394 ],  [ 179.78389441 ,  71.07953968], [ 179.94830596 ,  71.05551774]]

aCoords_gcs = np.array(aCoords_gcs)
aCoord_gcs_split = split_polygon_cross_idl(aCoords_gcs)