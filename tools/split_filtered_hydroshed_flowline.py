import os
from osgeo import ogr, osr, gdal
import numpy as np


from shapely.wkt import loads
from shapely.geometry import  MultiLineString, mapping, shape

from pyflowline.formats.convert_coordinates import convert_gcs_coordinates_to_flowline
from pyflowline.formats.export_flowline import export_flowline_to_geojson

aTopology=list()
aTopology_id=list()
aTopology_downid=list()
aFlowline = list()
aFlag_lake=list()
aLength=list()
lID_local = 0

sWorkspace_out  = '/compyfs/liao313/00raw/mesh/global/pyflowline'

def read_hydroshed_topology(sFilename_filtered_hydroshed_in):
    global aTopology_id
    global aTopology_downid
    global aFlowline 
    global aLength
    global aFlag_lake

    pDriver_shapefile = ogr.GetDriverByName('ESRI Shapefile')   
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
            return 1
        else:
            nUpstream = len(index[0])
            aUpstreamid = aTopology_id0[index[0]]
            aUpstreamindex = index[0]

            return nUpstream, aUpstreamid, aUpstreamindex
       
    def tag_upstream(lRiverID_in):
        global lID_local
        if(check_head_water(lRiverID_in)==1):            
            pass
        else:
            nUpstream, aUpstreamid, aUpstreamindex = find_upstream_flowline(lRiverID_in)
            if nUpstream > 0:                
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

    tag_upstream(lRiverID_in)
            

    return aFlowline_out

def split_filtered_hydroshed_flowline():
    global aFlowline
  
    global aTopology_id
    global aTopology_downid

    global aLength
    global aFlag_lake
    global lID_local

    sFilename_filtered_hydroshed = '/compyfs/liao313/00raw/mesh/global/filtered_12p5_2p5/filtered_12p5_2p5.shp'
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
        lNext_down = int(aTopology_downid[i])
        iFlag_lake = int(aFlag_lake[i])
        #other attributes
        dLength = float(aLength[i])        
        if lNext_down == 0: #this is a outlet, however, it may be associated with a lake instead of ocean
            if iFlag_lake == 0:   
                sBasin =    "{:05d}".format(iBasin)
                lID_local = 0
                aFlowline_basin =find_upstream(lRiverID)                                         
                #save as a geojson
                sFilename_json_in =  sWorkspace_out + '/basin' + sBasin+ '.geojson'
                export_flowline_to_geojson(aFlowline_basin, sFilename_json_in)
                iBasin = iBasin + 1                
            else:
                #this is a lake outlet
                pass
    

if __name__ == '__main__':
    split_filtered_hydroshed_flowline()