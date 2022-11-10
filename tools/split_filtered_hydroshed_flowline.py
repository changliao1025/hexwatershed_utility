import os
from osgeo import ogr, osr, gdal
import numpy as np


from shapely.wkt import loads
from shapely.geometry import  MultiLineString, mapping, shape

from pyflowline.formats.convert_coordinates import convert_gcs_coordinates_to_flowline

aTopology=list()
aTopology_id=list()
aTopology_downid=list()
aFlowline = list()
aFlag_lake=list()
aLength=list()


def read_hydroshed_topology(sFilename_filtered_hydroshed_in):
    global aTopology_id
    global aTopology_downid
    global aFlowline 
    global aLength
    global aFlag_lake

    pDriver_shapefile = ogr.GetDriverByName('ESRI Shapefile')   
    pDataset_shapefile = pDriver_shapefile.Open(sFilename_filtered_hydroshed_in, gdal.GA_ReadOnly)
    pLayer_shapefile = pDataset_shapefile.GetLayer(0)
    aFlowline_individual = list()
    
    nfeature_flowline = pLayer_shapefile.GetFeatureCount()  
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
            #for Line in aLine: 
            nLine = pGeometry_in.GetGeometryCount()
            for i in range(nLine):
                Line = pGeometry_in.GetGeometryRef(i)  
                dummy = loads( Line.ExportToWkt() )
                aCoords = dummy.coords                
                dummy1= np.array(aCoords)
                pLine = convert_gcs_coordinates_to_flowline(dummy1)
                pLine.lIndex = i
                aFlowline_individual.append(pLine)
                lID = lID 
        else:
            if sGeometry_type =='LINESTRING':
                dummy = loads( pGeometry_in.ExportToWkt() )
                aCoords = dummy.coords                
                dummy1= np.array(aCoords)
                pLine = convert_gcs_coordinates_to_flowline(dummy1)
                pLine.lIndex = 0
                aFlowline_individual.append(pLine)
                lID = lID 
            else:
                print(sGeometry_type)
                pass
        
        aFlowline.append(aFlowline_individual)        
        aTopology_id.append(lRiverID)
        aTopology_downid.append(lNext_down)
        aFlag_lake.append(iFlag_lake)
        aLength.append(dLength)

    #aTopology = np.array(aTopology)
    return aFlowline, aTopology_id, aTopology_downid

def find_upstream(lRiverID_in):
    global aTopology_id
    global aTopology_downid
    global aFlowline 
    global aLength
    global aFlag_lake

    dummy_index  = np.where(aTopology_downid == lRiverID_in )
    nUpstream = len(dummy_index)
    if nUpstream > 0:
        for j in range(nUpstream):
            
            lFlowlineID = aTopology_id[dummy_index[j]]
            

    global aFlowline

def split_filtered_hydroshed_flowline():
    global aFlowline
  
    global aTopology_id
    global aTopology_downid

    global aLength
    global aFlag_lake

    sFilename_filtered_hydroshed = '/compyfs/liao313/00raw/mesh/global/filtered_12p5_2p5/filtered_12p5_2p5.shp'
    sFilename_shapefile_in = sFilename_filtered_hydroshed

    if os.path.isfile(sFilename_shapefile_in):
        pass
    else:
        print('This shapefile does not exist: ', sFilename_shapefile_in )
        iReturn_code = 0
        return iReturn_code

    read_hydroshed_topology(sFilename_filtered_hydroshed)
  
   
    lID = 0

    #for pFeature_shapefile in pLayer_shapefile: 
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
                aFlowline_river = list()                             
                #then search in the lookup table
                iFlag_found_headwater = 0
                while iFlag_found_headwater  == 0 :
                    iFlag_found_headwater = find_upstream(lRiverID)      

                #save as a geojson

                
            else:
                #this is a lake outlet
                pass

    

if __name__ == '__main__':
    split_filtered_hydroshed_flowline()