import os, sys, stat
from pathlib import Path
import json
from netCDF4 import Dataset
from pyearth.system.define_global_variables import *    
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.ticker as mticker
import matplotlib.patches as mpatches
import cartopy.crs as ccrs
import matplotlib.cm as cm
import numpy as np
from shapely.wkt import loads
from osgeo import  osr, gdal, ogr
#desired_proj = ccrs.PlateCarree()
desired_proj = ccrs.Orthographic(central_longitude=-145, central_latitude=69, globe=None)
desired_proj = ccrs.Orthographic(central_longitude=  0.50*(-149.5+(-146.5)), \
        central_latitude= 0.50*(68.1+70.35), globe=None)

def map_mosart_variable(sFilename_out,sTitle, sUnit, aData, aExtent_in=None):

    sFilename_json0 = '/compyfs/liao313/04model/pyhexwatershed/sag/pyhexwatershed20220607001/pyflowline/mpas.geojson'
    
    sFilename_json = '/compyfs/liao313/04model/pyhexwatershed/sag/pyhexwatershed20220607001/hexwatershed/hexwatershed.json'
    
    if os.path.exists(sFilename_out):
        os.remove(sFilename_out)
    
    dData_min = 0.0
    dData_max = 1.0
    dData_max = np.percentile(aData, 95)
 

    fig = plt.figure( dpi=300 )
    fig.set_figwidth( 12 )
    fig.set_figheight( 12 )
    ax = fig.add_axes([0.1, 0.5, 0.75, 0.4] , projection=desired_proj )    
    ax.set_global()     
    
    dLat_min = 90
    dLat_max = -90
    dLon_min = 180
    dLon_max = -180
    pDriver = ogr.GetDriverByName('GeoJSON')
    pDataset = pDriver.Open(sFilename_json0, gdal.GA_ReadOnly)
    pLayer = pDataset.GetLayer(0)
    for pFeature in pLayer:
        pGeometry_in = pFeature.GetGeometryRef()
        sGeometry_type = pGeometry_in.GetGeometryName()
        lID =0 
        if sGeometry_type =='POLYGON':
            dummy0 = loads( pGeometry_in.ExportToWkt() )
            aCoords_gcs = dummy0.exterior.coords
            aCoords_gcs= np.array(aCoords_gcs)
            nvertex = len(aCoords_gcs)
            for i in range(nvertex):
                dLon = aCoords_gcs[i][0]
                dLat = aCoords_gcs[i][1]
                if dLon > dLon_max:
                    dLon_max = dLon
                if dLon < dLon_min:
                    dLon_min = dLon
                if dLat > dLat_max:
                    dLat_max = dLat
                if dLat < dLat_min:
                    dLat_min = dLat
            polygon = mpatches.Polygon(aCoords_gcs[:,0:2], closed=True, linewidth=0.25, \
                alpha=0.8, edgecolor = 'black',facecolor='none', \
                    transform=ccrs.PlateCarree() )
            ax.add_patch(polygon)  

    cmap = cm.get_cmap('Spectral')
    cmap_reversed = cmap.reversed()
    norm=plt.Normalize(dData_min,dData_max)
    aData[np.where(aData > dData_max)] = dData_max
    with open(sFilename_json) as json_file:
        data = json.load(json_file)     
        ncell = len(data)
        lID =0 
        for i in range(ncell):
            pcell = data[i]
            lCellID = int(pcell['lCellID'])
            lCellID_downslope = int(pcell['lCellID_downslope'])
            x_start=float(pcell['dLongitude_center_degree'])
            y_start=float(pcell['dLatitude_center_degree'])
            dfac = float(pcell['DrainageArea'])
            dummy = aData[i]
            avertex = pcell['vVertex']
            nvertex = len(avertex)
            aLocation= np.full( (nvertex, 2), 0.0, dtype=float )
            #this is the cell
            #get the vertex
            for k in range(nvertex):
                aLocation[k,0] = avertex[k]['dLongitude_degree']
                aLocation[k,1] = avertex[k]['dLatitude_degree']
                #if aLocation[k,0] > dLon_max:
                #    dLon_max = aLocation[k,0]
                #if aLocation[k,0] < dLon_min:
                #    dLon_min = aLocation[k,0]
                #if aLocation[k,1] > dLat_max:
                #    dLat_max = aLocation[k,1]
                #if aLocation[k,1] < dLat_min:
                #    dLat_min = aLocation[k,1]
            color_index = (dummy-dData_min ) /(dData_max - dData_min )
            rgba = cmap_reversed(color_index)
            polygon = mpatches.Polygon(aLocation, closed=True, linewidth=0.3,\
                facecolor=rgba,\
                edgecolor='none',transform=ccrs.PlateCarree() )
            
            ax.add_patch(polygon)     
    #trasform elevation
    sm = plt.cm.ScalarMappable(cmap=cmap_reversed, norm=norm)
    
    sm.set_array(aData)
    cb = fig.colorbar(sm, ax=ax)        
    
    cb.ax.get_yaxis().set_ticks_position('right')
    cb.ax.get_yaxis().labelpad = 5
    cb.ax.set_ylabel(sUnit, rotation=90)
    cb.ax.get_yaxis().set_label_position('left')
    cb.ax.tick_params(labelsize=6) 
    
    
    if aExtent_in is not None:
        aExtent = aExtent_in #= [-78.5,-75.5, 39.2,42.5]
    else:
        aExtent = [dLon_min, dLon_max, dLat_min, dLat_max]
    ax.set_extent(aExtent)
    ax.coastlines()#resolution='110m')        
    #crs=ccrs.PlateCarree(),
    gl = ax.gridlines( draw_labels=True,\
                  linewidth=0.2, color='gray', alpha=0.3, linestyle='--')
    gl.xlabel_style = {'size': 8, 'color': 'k', 'rotation':0, 'ha':'right'}
    gl.ylabel_style = {'size': 8, 'color': 'k', 'rotation':90,'weight':'normal'}
    ax.set_title(sTitle , loc='center')
    
    sText = ' ' 
    
    ax.text(0.05, 0.95, sText, \
    verticalalignment='top', horizontalalignment='left',\
            transform=ax.transAxes, \
            color='black', fontsize=8)
    sText = 'Mesh type: ' + 'MPAS' 
    ax.text(0.05, 0.90, sText, \
    verticalalignment='top', horizontalalignment='left',\
            transform=ax.transAxes, \
            color='black', fontsize=8)
    
    plt.savefig(sFilename_out, bbox_inches='tight')
    plt.close(fig)
    
    
def map_mosart_output():

    iYear_start = 1979
    iYear_end = 1980
    iMonth_start=1
    iMonth_end=12

    sWorkspace_simulation_case_run = '/compyfs/liao313/e3sm_scratch/e3sm20220701048/run'
    sCase = 'e3sm20220701048'
    
    sWorkspace_analysis_case = '/compyfs/liao313/04model/e3sm/sag/analysis/e3sm20220701048/'

    sVariable_discharge  = 'RIVER_DISCHARGE_OVER_LAND_LIQ'
    sTitle='River discharge'
    sUnit = r'Units: $m^{3} s^{-1}$'
    if not os.path.exists(sWorkspace_analysis_case):
        os.makedirs(sWorkspace_analysis_case)  

    aExtent_in = None #[-153,-143, 65,75]

    for iYear in range(iYear_start, iYear_end+1):

        sYear =   "{:04d}".format(iYear)
        for iMonth in range(iMonth_start, iMonth_end + 1):
            sMonth = str(iMonth).zfill(2)           
    
            sDummy = '.mosart.h0.' + sYear + '-' + sMonth + sExtension_netcdf
            sFilename = sWorkspace_simulation_case_run + slash + sCase + sDummy
            print(sFilename)
            if os.path.exists(sFilename):
                #print("Yep, I can read that file: " + sFilename)                
                pass
            else:
                print(sFilename + ' is missing')
                print("Nope, the path doesn't reach your file. Go research filepath in python")
                return
    
            aDatasets = Dataset(sFilename)
            sFilename_out = sWorkspace_analysis_case + slash + 'discharge' + sYear + sMonth + '.png'
            
            for sKey, aValue in aDatasets.variables.items():
                if sKey.lower() == sVariable_discharge.lower() :

                    aData_out = ((aValue[:]).data)[0]     

                    map_mosart_variable(sFilename_out,sTitle, sUnit, aData_out, aExtent_in=aExtent_in)
                    break
                    

        

    

    return

if __name__ == '__main__':
    map_mosart_output()