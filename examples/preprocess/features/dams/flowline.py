if self.iFlag_dam == 1:
            sFilename_flowline_filter = self.sFilename_flowline_filter_geojson
            aFlowline_basin_filtered_raw, pProjection_geojson = read_flowline_geojson( sFilename_flowline_filter )
            aVertex_filtered = find_flowline_vertex(aFlowline_basin_filtered_raw)

            ptimer.start()
            nFlowline_before = len(aFlowline_basin_filtered_raw)
            sFilename_dam = self.sFilename_dam
            #obtain dam lookup table C
            aData_dam = text_reader_string(sFilename_dam, iSkipline_in =1,cDelimiter_in=',' )
            sFilename_flowline_topo = self.sFilename_flowline_topo
            #obtain whole topology B
            aData_flowline_topo = text_reader_string(sFilename_flowline_topo, iSkipline_in =1,cDelimiter_in=',' )
            aFromFlowline = aData_flowline_topo[:,1].astype(int).ravel()
            aToFlowline = aData_flowline_topo[:,2].astype(int).ravel()
            sFilename_flowline_raw = self.sFilename_flowline_raw
            #find A lookup table
            aNHDPlusID_filter = read_nhdplus_flowline_geojson_attribute(sFilename_flowline_filter)
            ndam = len(aData_dam)
            aNHDPlusID_dams_headwater = list()
            aFlowline_dams_nonheadwater = list()
            aVertex_dams_nonheadwater = list()
            n=0
            for j in range(0, ndam):
                dLon = float(aData_dam[j][1])
                dLat = float(aData_dam[j][0])
                sDam = aData_dam[j][4]
                #individual ID
                lNHDPlusID = int(aData_dam[j][5])
                #if lNHDPlusID in aNHDPlusID_filter:
                if lNHDPlusID in aNHDPlusID_filter: #this is already included in A
                    #change flag by id
                    for k in range(len(aFlowline_basin_filtered_raw)):
                        #remove this flowline by ID
                        if aFlowline_basin_filtered_raw[k].lNHDPlusID == lNHDPlusID:
                            aFlowline_basin_filtered_raw[k].iFlag_dam =1
                            break
                    pass
                else:
                    aNHDPlusID_dams_headwater.append(lNHDPlusID)
                    #not in A, so we need to trace it down

                    aNHDPlusID_dam_nonheadwater = track_nhdplus_flowline(aNHDPlusID_filter, aFromFlowline, aToFlowline, lNHDPlusID)
                    aFlowline_dam_nonheadwater = extract_nhdplus_flowline_shapefile_by_attribute(sFilename_flowline_raw, aNHDPlusID_dam_nonheadwater )
                    #clean up
                    aVertex_dam_nonheadwater = find_flowline_vertex(aFlowline_dam_nonheadwater)
                    dThreshold = 1.0
                    iFlag_found=0
                    n=n+len(aFlowline_dam_nonheadwater)
                    #print(j,n)
                    for pVertex in aVertex_filtered:
                        if iFlag_found ==1:
                            break

                        for i in range( len(aVertex_dam_nonheadwater) ):
                            pVertex_dam = aVertex_dam_nonheadwater[i]
                            dDistance = pVertex.calculate_distance(pVertex_dam)
                            if  dDistance<= dThreshold:
                                aVertex_dam_nonheadwater[i] = pVertex
                                for k in range(len(aFlowline_dam_nonheadwater)):
                                    if aFlowline_dam_nonheadwater[k].pVertex_end == pVertex_dam:
                                        #update
                                        aEdge = aFlowline_dam_nonheadwater[k].aEdge
                                        aEdge[-1] = pyedge( aEdge[-1].pVertex_start, pVertex)
                                        aFlowline_dam_nonheadwater[k] = pyflowline(aEdge)

                                iFlag_found = 1
                                break
                            else:
                                #print(dDistance)
                                pass

                    aVertex_dams_nonheadwater.append(aVertex_dam_nonheadwater)

                    aFlowline_dams_nonheadwater.append(aFlowline_dam_nonheadwater)


            aFlowline_dams_headwater = extract_nhdplus_flowline_shapefile_by_attribute(sFilename_flowline_raw, aNHDPlusID_dams_headwater )
            for i in range(len(aFlowline_dams_headwater)):
                aFlowline_dams_headwater[i].iFlag_dam = 1

            aFlowline_dams_nonheadwater_all = [item for sublist in aFlowline_dams_nonheadwater for item in sublist]
            aVertex_dam_nonheadwater_all = [item for sublist in aVertex_dams_nonheadwater for item in sublist]

            aFlowline_basin_filtered = aFlowline_basin_filtered_raw + aFlowline_dams_headwater + aFlowline_dams_nonheadwater_all
            aVertex_dam = find_flowline_vertex(aFlowline_dams_headwater)

            if self.iFlag_debug ==1:
                sFilename_out = 'flowline_vertex_filtered.geojson'
                sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
                export_vertex_to_geojson( aVertex_filtered, sFilename_out)
                sFilename_out = 'flowline_vertex_dam.geojson'
                sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
                export_vertex_to_geojson( aVertex_dam, sFilename_out)
                sFilename_out = 'flowline_vertex_dam_nonheadwater.geojson'
                sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
                export_vertex_to_geojson( aVertex_dam_nonheadwater_all, sFilename_out)

            nFlowline_after = len(aFlowline_basin_filtered)
            print('Basin ',  self.sBasinID, ' has dam', nFlowline_before, nFlowline_after)
            ptimer.stop()
        else:
            print('Basin ',  self.sBasinID, ' has no dam')
            sFilename_flowline_filter = self.sFilename_flowline_filter_geojson #sFilename_flowline_filter = self.sFilename_flowline_filter
            aFlowline_basin_filtered, pProjection_geojson = read_flowline_geojson( sFilename_flowline_filter )
            #aVertex_filtered = find_flowline_vertex(aFlowline_basin_filtered)
            pass