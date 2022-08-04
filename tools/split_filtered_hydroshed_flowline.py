

def split_filtered_hydroshed_flowline():

    sFilename_filtered_hydroshed = '/compyfs/liao313/00raw/mesh/global/filtered_12p5_2p5/filtered_12p5_2p5.shp'
    sFilename_shapefile_in = sFilename_filtered_hydroshed

    if os.path.isfile(sFilename_shapefile_in):
        pass
    else:
        print('This shapefile does not exist: ', sFilename_shapefile_in )
        iReturn_code = 0
        return iReturn_code

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

    for pFeature_shapefile in pLayer_shapefile:        
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
                    
                    for Line in aLine: 
                        dummy = loads( Line.ExportToWkt() )
                        aCoords = dummy.coords                
                        dummy1= np.array(aCoords)
                        pLine = convert_gcs_coordinates_to_flowline(dummy1)
                        pLine.lIndex = lID
                        pLine.lNHDPlusID= lNHDPlusID
                        aFlowline.append(pLine)
                        lID = lID + 1

                else:
                    if sGeometry_type =='LINESTRING':
                        dummy = loads( pGeometry_in.ExportToWkt() )
                        aCoords = dummy.coords                
                        dummy1= np.array(aCoords)
                        pLine = convert_gcs_coordinates_to_flowline(dummy1)
                        pLine.lIndex = lID
                        pLine.lNHDPlusID= lNHDPlusID
                        aFlowline.append(pLine)
                        lID = lID + 1

                    else:
                        print(sGeometry_type)
                        pass

            else:
                #this is a lake outlet

        

if __name__ == '__main__':
    split_filtered_hydroshed_flowline()