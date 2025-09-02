import os
import numpy as np
from osgeo import gdal, osr, ogr

def map_wkt_polygon(pPolygon_wkt, sFilename_png_out):
    #this is optional, if you want to plot the polygon, you can use this function, it calls some functions from the pyearth package
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon
    from matplotlib.collections import PatchCollection
    import matplotlib as mpl
    import cartopy.crs as ccrs
    from cartopy.io.img_tiles import OSM
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
    from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates
    from pyearth.visual.map.map_servers import calculate_zoom_level, calculate_scale_denominator
    from pyearth.visual.map.zebra_frame import zebra_frame
    pSRS_wgs84 = ccrs.PlateCarree()  # for latlon data only
    #get the list of coordinates in the wktpolygon
    pPolygon = ogr.CreateGeometryFromWkt(pPolygon_wkt)
    aCoords_gcs = get_geometry_coordinates(pPolygon)
    aLon = aCoords_gcs[:,0]
    aLat = aCoords_gcs[:,1]
    dLon_max = np.max(aLon)
    dLon_min = np.min(aLon)
    dLat_max = np.max(aLat)
    dLat_min = np.min(aLat)
    pProjection_map = ccrs.Orthographic(central_longitude=0.50*(
            dLon_max+dLon_min),  central_latitude=0.50*(dLat_max+dLat_min), globe=None)

    fig = plt.figure(dpi=150)
    fig.set_figwidth(8)
    fig.set_figheight(8)
    ax = fig.add_axes([0.08, 0.1, 0.62, 0.7], projection=pProjection_map)
    image_size = [1000, 1000]
    aExtent = [dLon_min, dLon_max, dLat_min, dLat_max]
    scale_denominator = calculate_scale_denominator(aExtent, image_size)
    pSrc = osr.SpatialReference()
    pSrc.ImportFromEPSG(3857) # mercator
    pProjection = pSrc.ExportToWkt()
    iBasemap_zoom_level = calculate_zoom_level(scale_denominator, pProjection)
    osm_tiles = OSM()
    #Add the OSM image to the map
    ax.add_image(osm_tiles, iBasemap_zoom_level)
    aPolygon = list()
    aPolygon.append(aCoords_gcs)
    aPatch = [Polygon(poly, closed=True) for poly in aPolygon]
    pPC = PatchCollection(aPatch, alpha=0.5,
                                  edgecolor='black',
                                  facecolor='none',
                                  linewidths=2,
                                  transform=pSRS_wgs84)
    ax.add_collection(pPC)
    marginx = (dLon_max - dLon_min) / 20
    marginy = (dLat_max - dLat_min) / 20
    if (dLat_max + marginy)> 90:
        dLat_max = 90
    else:
        dLat_max = dLat_max + marginy
    if (dLat_min - marginy) < -90:
        dLat_min = -90
    else:
        dLat_min = dLat_min - marginy
    if (dLon_max + marginx) > 180:
        dLon_max = 180
    else:
        dLon_max = dLon_max + marginx
    if (dLon_min - marginx) < -180:
        dLon_min = -180
    else:
        dLon_min = dLon_min - marginx
    aExtent_extend = [dLon_min, dLon_max, dLat_min, dLat_max]
    minx, maxx, miny, maxy = aExtent_extend
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=1, color='gray', alpha=0.5, linestyle='--',
                      xlocs=np.arange(minx, maxx+(maxx-minx)/9, (maxx-minx)/8),
                      ylocs=np.arange(miny, maxy+(maxy-miny)/9, (maxy-miny)/8))
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlocator = mpl.ticker.MaxNLocator(4)
    gl.ylocator = mpl.ticker.MaxNLocator(4)
    gl.xlabel_style = {'size': 10, 'color': 'k', 'rotation': 0, 'ha': 'right'}
    gl.ylabel_style = {'size': 10, 'color': 'k',
                       'rotation': 90, 'weight': 'normal'}

    ax.set_extent(aExtent_extend, crs = pSRS_wgs84)
    ax.zebra_frame(crs=pSRS_wgs84, iFlag_outer_frame_in=1)
    plt.savefig(sFilename_png_out, bbox_inches='tight')
    print('The plot is saved to: ', sFilename_png_out)

    #save the polygon as a geojson file
    sFilename_geojson = sFilename_png_out.replace('.png', '.geojson')
    pDriver = ogr.GetDriverByName('GeoJSON')
    if os.path.exists(sFilename_geojson):
        pDriver.DeleteDataSource(sFilename_geojson)

    pDS = pDriver.CreateDataSource(sFilename_geojson)

    pLayer = pDS.CreateLayer('polygon', srs=None)
    pFieldDefn = ogr.FieldDefn('id', ogr.OFTInteger)
    pLayer.CreateField(pFieldDefn)
    pFeature = ogr.Feature(pLayer.GetLayerDefn())
    pFeature.SetField('id', 1)
    pFeature.SetGeometry(pPolygon)
    pLayer.CreateFeature(pFeature)
    pFeature = None
    pDS = None
    print('The polygon is saved to: ', sFilename_geojson)
    #save each point
    sFilename_geojson = sFilename_png_out.replace('.png', '_point.geojson')
    if os.path.exists(sFilename_geojson):
        pDriver.DeleteDataSource(sFilename_geojson)

    pDS = pDriver.CreateDataSource(sFilename_geojson)
    pLayer = pDS.CreateLayer('point', srs=None)
    pFieldDefn = ogr.FieldDefn('id', ogr.OFTInteger)
    pLayer.CreateField(pFieldDefn)
    #add lon and lat to the point
    pFieldDefn = ogr.FieldDefn('lon', ogr.OFTReal)
    pLayer.CreateField(pFieldDefn)
    pFieldDefn = ogr.FieldDefn('lat', ogr.OFTReal)
    pLayer.CreateField(pFieldDefn)
    npoint = len(aCoords_gcs)
    for i in range(npoint-1):
        pPoint = ogr.Geometry(ogr.wkbPoint)
        pPoint.AddPoint(aCoords_gcs[i,0], aCoords_gcs[i,1])
        pFeature = ogr.Feature(pLayer.GetLayerDefn())
        pFeature.SetField('id', i)
        pFeature.SetField('lon', aCoords_gcs[i,0])
        pFeature.SetField('lat', aCoords_gcs[i,1])
        pFeature.SetGeometry(pPoint)
        pLayer.CreateFeature(pFeature)
        pFeature = None
        pass
    pDS = None
    print('The points are saved to: ', sFilename_geojson)

    return