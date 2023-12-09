import os, sys, stat
from osgeo import ogr, osr, gdal
import numpy as np

from pathlib import Path
from os.path import realpath


from shapely.wkt import loads


from pyflowline.formats.convert_coordinates import convert_gcs_coordinates_to_flowline
from pyflowline.formats.export_flowline import export_flowline_to_geojson


sys.setrecursionlimit(10000)

aTopology=list()
aTopology_id=list()
aTopology_downid=list()
aFlowline = list()
aFlag_lake=list()
aLength=list()
lID_local = 0

#===========================
#setup workspace path
#===========================
sPath_parent = str(Path(__file__).parents[2]) # data is located two dir's up
sPath_data = realpath( sPath_parent +  '/data/' )
#sWorkspace_input =  str(Path(sPath_data)  /  'input')
sWorkspace_output = sPath_parent +  '/data/asia/pyflowline'


class TailRecurseException(Exception):
    def __init__(self, args, kwargs):
        self.args = args
        self.kwargs = kwargs

def tail_call_optimized(g):
    """
    This function decorates a function with tail call
    optimization. It does this by throwing an exception
    if it is it's own grandparent, and catching such
    exceptions to fake the tail call optimization.
    
    This function fails if the decorated
    function recurses in a non-tail context.
    """
    def func(*args, **kwargs):
        f = sys._getframe()
        if f.f_back and f.f_back.f_back  and f.f_back.f_back.f_code == f.f_code:
            raise TailRecurseException(args, kwargs)
        else:
            while 1:
                try:
                    return g(*args, **kwargs)
                except TailRecurseException as e:
                    args = e.args
                    kwargs = e.kwargs

    func.__doc__ = g.__doc__
    return func

def read_hydroshed_topology(sFilename_filtered_hydroshed_in):
    global aTopology_id
    global aTopology_downid
    global aFlowline 
    global aLength
    global aFlag_lake

    #pDriver_shapefile = ogr.GetDriverByName('ESRI Shapefile')   
    pDriver_shapefile = ogr.GetDriverByName('GeoJSON')  
    pDataset_shapefile = pDriver_shapefile.Open(sFilename_filtered_hydroshed_in, gdal.GA_ReadOnly)
    pLayer_shapefile = pDataset_shapefile.GetLayer(0)
    
    
    nfeature_flowline = pLayer_shapefile.GetFeatureCount()  
    lID = 0
    for i in range(nfeature_flowline):   
        pFeature_shapefile = pLayer_shapefile.GetFeature(i)
        pGeometry_in = pFeature_shapefile.GetGeometryRef()
        sGeometry_type = pGeometry_in.GetGeometryName()

        #get the flowline attribute, we only accept a flowline if it flows into the ocean
        lRiverID = int(pFeature_shapefile.GetField("HYRIV_ID"))        
        lNext_down = int(pFeature_shapefile.GetField("NEXT_DOWN"))        
        iFlag_lake = int(pFeature_shapefile.GetField("ENDORHEIC"))
        dLength = float(pFeature_shapefile.GetField("LENGTH_KM"))
        if(sGeometry_type == 'MULTILINESTRING'):                    
            aFlowline_individual = list()
            #for Line in aLine: 
            nLine = pGeometry_in.GetGeometryCount()
            for i in range(nLine):
                Line = pGeometry_in.GetGeometryRef(i)  
                dummy = loads( Line.ExportToWkt() )
                aCoords = dummy.coords                
                dummy1= np.array(aCoords)
                pLine0 = convert_gcs_coordinates_to_flowline(dummy1)                
                aFlowline_individual.append(pLine0)
                pass
            #merge?
            pFlowline = pLine0
            pFlowline.lIndex = lID
            pFlowline.lFlowlineID = lID     
        else:
            if sGeometry_type =='LINESTRING':
                dummy = loads( pGeometry_in.ExportToWkt() )
                aCoords = dummy.coords                
                dummy1= np.array(aCoords)
                pFlowline = convert_gcs_coordinates_to_flowline(dummy1)
                pFlowline.lIndex = lID
                pFlowline.lFlowlineID = lID
                
                
            else:
                print(sGeometry_type)
                pass
        lID = lID + 1        
        aFlowline.append(pFlowline)        
        aTopology_id.append(lRiverID)
        aTopology_downid.append(lNext_down)
        aFlag_lake.append(iFlag_lake)
        aLength.append(dLength)

    #aTopology = np.array(aTopology)
    #sFilename_json_in =  '/compyfs/liao313/00raw/mesh/global/flowline.geojson'
    #export_flowline_to_geojson(aFlowline, sFilename_json_in)
    return aFlowline, aTopology_id, aTopology_downid

def find_upstream(lRiverID_in):
    global aTopology_id
    global aTopology_downid
    global aFlowline 
    global aLength
    global aFlag_lake

    global lID_local

    aFlowline_out=list()
    def check_head_water(lRiverID_in):
        aTopology_downid0 = np.array(aTopology_downid).ravel()
        index = np.where(aTopology_downid0 == lRiverID_in)
        if len(index[0]) ==0 :
            return 1
        else:
            return 0

    def find_upstream_flowline(lRiverID_in):
        aTopology_id0 = np.array(aTopology_id).ravel()
        aTopology_downid0 = np.array(aTopology_downid).ravel()
        index = np.where(aTopology_downid0 == lRiverID_in)
        if len(index[0]) ==0 :            
            nUpstream=0
            aUpstreamid= None
            aUpstreamindex = None
        else:
            nUpstream = len(index[0])
            aUpstreamid = aTopology_id0[index[0]]
            aUpstreamindex = index[0]

        return nUpstream, aUpstreamid, aUpstreamindex
    
    #@tail_call_optimized
    def tag_upstream(lRiverID_in):
        global lID_local
        if(check_head_water(lRiverID_in)==1):            
            pass
        else:
            nUpstream, aUpstreamid, aUpstreamindex = find_upstream_flowline(lRiverID_in)
            if nUpstream > 0:                
                if nUpstream == 1:
                    pFlowline = aFlowline[ aUpstreamindex[0] ]                    
                    pFlowline.lFlowlineID = lID_local
                    aFlowline_out.append(pFlowline)
                    lID_local = lID_local + 1
                    tag_upstream(  aUpstreamid[0]  )          
                else:
                    for j in range(nUpstream):
                        pFlowline = aFlowline[ aUpstreamindex[j] ]                    
                        pFlowline.lFlowlineID = lID_local
                        aFlowline_out.append(pFlowline)
                        lID_local = lID_local + 1
                        tag_upstream(  aUpstreamid[j]  )            

                    pass
            else:
                pass
    
    lID_local = 0

    index = np.where(np.array(aTopology_id).ravel() == lRiverID_in)
    pFlowline = aFlowline[ index[0][0] ]
    aFlowline_out.append(pFlowline)

    tag_upstream(lRiverID_in)            

    return aFlowline_out

def split_filtered_hydroshed_flowline():
    global aFlowline
  
    global aTopology_id
    global aTopology_downid

    global aLength
    global aFlag_lake
    global lID_local

    sFilename_filtered_hydroshed = sPath_parent  + '/data/conus/hydroshed_conus.geojson'
    sFilename_filtered_hydroshed = '/compyfs/liao313/00raw/mesh/southamerica/southamerica.geojson'
    sFilename_filtered_hydroshed = '/compyfs/liao313/00raw/mesh/asia/china.geojson'
    if os.path.isfile(sFilename_filtered_hydroshed):
        pass
    else:
        print('This shapefile does not exist: ', sFilename_filtered_hydroshed )
        iReturn_code = 0
        return iReturn_code   

    read_hydroshed_topology(sFilename_filtered_hydroshed)  
   
    lID = 0
    iBasin = 1    
    nfeature_flowline = len(aFlowline)      
    for i in range(nfeature_flowline):  
        #get the flowline attribute, we only accept a flowline if it flows into the ocean
        lRiverID = int(aTopology_id[i])   
        #40613666 yangzte river
        #60443230 amazon
        if lRiverID == 40613666: #
            print('l')
        else:
            continue
        lNext_down = int(aTopology_downid[i])
        iFlag_lake = int(aFlag_lake[i])
        #other attributes
        dLength = float(aLength[i])        
        if lNext_down == 0: #this is a outlet, however, it may be associated with a lake instead of ocean
            if iFlag_lake == 0:   
                sBasin =    "{:05d}".format(iBasin)
                sRiver =  "{:0d}".format(lRiverID)
                print(sRiver)
                lID_local = 0
                aFlowline_basin =find_upstream(lRiverID)                                         
                #save as a geojson
                sFilename_json_in =  sWorkspace_output + '/river' + sRiver+ '.geojson'
                export_flowline_to_geojson(aFlowline_basin, sFilename_json_in)
                iBasin = iBasin + 1                
            else:
                #this is a lake outlet
                pass    

if __name__ == '__main__':
    split_filtered_hydroshed_flowline()