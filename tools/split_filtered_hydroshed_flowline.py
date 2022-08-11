import os
from osgeo import ogr, osr, gdal
import numpy as np


from shapely.wkt import loads
from shapely.geometry import  MultiLineString, mapping, shape

from pyflowline.formats.convert_coordinates import convert_gcs_coordinates_to_flowline



def read_hydroshed_topology(sFilename_filtered_hydroshed_in):
    aTopology= list()
    pDriver_shapefile = ogr.GetDriverByName('ESRI Shapefile')   
    pDataset_shapefile = pDriver_shapefile.Open(sFilename_filtered_hydroshed_in, gdal.GA_ReadOnly)
    pLayer_shapefile = pDataset_shapefile.GetLayer(0)
    
    
    nfeature_flowline = pLayer_shapefile.GetFeatureCount()  
    for i in range(nfeature_flowline):   
        pFeature_shapefile = pLayer_shapefile.GetFeature(i)     
        pGeometry_in = pFeature_shapefile.GetGeometryRef()
        sGeometry_type = pGeometry_in.GetGeometryName()

        #get the flowline attribute, we only accept a flowline if it flows into the ocean
        lRiverID = int(pFeature_shapefile.GetField("HYRIV_ID"))        
        lNext_down = int(pFeature_shapefile.GetField("NEXT_DOWN"))

        pTopology = np.array([lRiverID,lNext_down])
        aTopology.append(pTopology)

    aTopology = np.array(aTopology)
    return aTopology

def split_filtered_hydroshed_flowline():

    sFilename_filtered_hydroshed = '/compyfs/liao313/00raw/mesh/global/filtered_12p5_2p5/filtered_12p5_2p5.shp'
    sFilename_shapefile_in = sFilename_filtered_hydroshed

    if os.path.isfile(sFilename_shapefile_in):
        pass
    else:
        print('This shapefile does not exist: ', sFilename_shapefile_in )
        iReturn_code = 0
        return iReturn_code

    aTopology = read_hydroshed_topology(sFilename_filtered_hydroshed)
    aTopology_id = aTopology[:, 0].flatten()
    aTopology_downid = aTopology[:, 1].flatten()
    

    aFlowline=list()
    pDriver_shapefile = ogr.GetDriverByName('ESRI Shapefile')   
    pDataset_shapefile = pDriver_shapefile.Open(sFilename_shapefile_in, gdal.GA_ReadOnly)
    pLayer_shapefile = pDataset_shapefile.GetLayer(0)
    pSpatialRef_shapefile = pLayer_shapefile.GetSpatialRef()
    pSpatial_reference_gcs = osr.SpatialReference()
    pSpatial_reference_gcs.ImportFromEPSG(4326)
    pSpatial_reference_gcs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    comparison = pSpatialRef_shapefile.IsSame(pSpatial_reference_gcs)
    if(comparison != 1):
        iFlag_transform =1
        pTransform = osr.CoordinateTransformation(pSpatialRef_shapefile, pSpatial_reference_gcs)
    else:
        iFlag_transform =0
   
    lID = 0

    #for pFeature_shapefile in pLayer_shapefile: 
    nfeature_flowline = pLayer_shapefile.GetFeatureCount()       
    for i in range(nfeature_flowline):   
        pFeature_shapefile = pLayer_shapefile.GetFeature(i)
        pGeometry_in = pFeature_shapefile.GetGeometryRef()
        sGeometry_type = pGeometry_in.GetGeometryName()

        #get the flowline attribute, we only accept a flowline if it flows into the ocean
        lRiverID = int(pFeature_shapefile.GetField("HYRIV_ID"))        
        lNext_down = int(pFeature_shapefile.GetField("NEXT_DOWN"))
        iFlag_lake = int(pFeature_shapefile.GetField("ENDORHEIC"))
        #other attributes
        dLength = float(pFeature_shapefile.GetField("LENGTH_KM"))

        if (iFlag_transform ==1): #projections are different
            pGeometry_in.Transform(pTransform)
        if (pGeometry_in.IsValid()):
            pass
        else:
            print('Geometry issue')

        if lNext_down == 0: #this is a outlet, however, it may be associated with a lake instead of ocean

            if iFlag_lake == 0: 
                #this is an ocean outlet
                if(sGeometry_type == 'MULTILINESTRING'):
                    nLine = pGeometry_in.GetGeometryCount()
                    for i in range(nLine):
                    #for Line in aLine:
                        Line = pGeometry_in.GetGeometryRef(i)  
                        dummy = loads( Line.ExportToWkt() )
                        aCoords = dummy.coords                
                        dummy1= np.array(aCoords)
                        pLine = convert_gcs_coordinates_to_flowline(dummy1)
                        pLine.lIndex = lID
                 
                        aFlowline.append(pLine)
                        lID = lID + 1

                else:
                    if sGeometry_type =='LINESTRING':
                        dummy = loads( pGeometry_in.ExportToWkt() )
                        aCoords = dummy.coords                
                        dummy1= np.array(aCoords)
                        pLine = convert_gcs_coordinates_to_flowline(dummy1)
                        pLine.lIndex = lID
                        
                        aFlowline.append(pLine)
                        lID = lID + 1

                    else:
                        print(sGeometry_type)
                        pass
                
                #then search in the lookup table
                iFlag_found_headwater = 0
                while iFlag_found_headwater  == 0 :
                    dummy_index  = np.where(aTopology_downid == lRiverID )
                    nUpstream = len(dummy_index)
                    if nUpstream > 0:
                        for j in range(nUpstream):
                            
                        lFlowlineID = aTopology_id[dummy_index[0]]
                        pFeature_shapefile = pLayer_shapefile.GetFeature(dummy_index[0])
                        pGeometry_in = pFeature_shapefile.GetGeometryRef()
                        sGeometry_type = pGeometry_in.GetGeometryName()
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
                                pLine.lIndex = lID
                                
                                aFlowline.append(pLine)
                                lID = lID + 1

                        else:
                            if sGeometry_type =='LINESTRING':
                                dummy = loads( pGeometry_in.ExportToWkt() )
                                aCoords = dummy.coords                
                                dummy1= np.array(aCoords)
                                pLine = convert_gcs_coordinates_to_flowline(dummy1)
                                pLine.lIndex = lID
                               
                                aFlowline.append(pLine)
                                lID = lID + 1

                            else:
                                print(sGeometry_type)
                                pass

                    else:
                        iFlag_found_headwater = 1




            else:
                #this is a lake outlet
                pass



    return aFlowline     

if __name__ == '__main__':
    split_filtered_hydroshed_flowline()