import os
import math
import os
import datetime
import textwrap
import numpy as np
from osgeo import osr
from urllib.error import URLError
import matplotlib as mpl
import matplotlib.pyplot as plt

from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon as mpolygon
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import shapely.geometry as sgeom
from pyearth.system.define_global_variables import sPath_proj_lib
from pyearthviz.map.zebra_frame import zebra_frame
from pyearth.gis.geometry.convert_longitude_range import convert_360_to_180
from pyearthviz.formatter import OOMFormatter
from pyearthviz.map import RasterTileServer

iYear_current = datetime.datetime.now().year
sYear = str(iYear_current)

def map_mpas_mesh_file(sFilename_mpas_mesh_netcdf_in,
    sFilename_output_in=None,
    iFlag_color_in=None,
    iFlag_zebra_in=None,
    iFlag_fill_in=None,
    iFont_size_in=None,
    aBasemap_provider_in=None,
    iBasemap_zoom_level_in=None,
    sColor_in=None,
    sTitle_in=None,
    iDPI_in=None,
    iSize_x_in=None,
    iSize_y_in=None,
    dLongitude_center_in=None,
    dLatitude_center_in=None,
    dLinewidth_in=None,
    sFont_in=None,
    aLegend_in=None,
    aExtent_in=None,
    pProjection_map_in=None,
    pProjection_data_in=None,):
    import netCDF4 as nc
    pSRS_wgs84 = ccrs.PlateCarree()  # for latlon data only
    pSRS_geodetic = ccrs.Geodetic()

    if os.path.exists(sFilename_mpas_mesh_netcdf_in) is False:
        print("File does not exist")
        return

    if iDPI_in is not None:
        iDPI = iDPI_in
    else:
        iDPI = 150

    if iSize_x_in is not None:
        iSize_x = iSize_x_in
    else:
        iSize_x = 8

    if iSize_y_in is not None:
        iSize_y = iSize_y_in
    else:
        iSize_y = 8

    if iFont_size_in is not None:
        iFont_size = iFont_size_in
    else:
        iFont_size = 12

    if iFlag_color_in is not None:
        iFlag_color = iFlag_color_in
    else:
        iFlag_color = 0

    if iFlag_fill_in is not None:
        iFlag_fill = iFlag_fill_in
    else:
        iFlag_fill = False

    if iFlag_zebra_in is not None:
        iFlag_zebra = iFlag_zebra_in
    else:
        iFlag_zebra = 0

    if dLongitude_center_in is not None:
        dLongitude_center = dLongitude_center_in
    else:
        dLongitude_center = 0.0

    if dLatitude_center_in is not None:
        dLatitude_center = dLatitude_center_in
    else:
        dLatitude_center = 0.0

    if aExtent_in is not None:
        aExtent = aExtent_in
        minx, maxx,miny,  maxy = aExtent
    else:
        #because the earth is roung, we can set it as the half of the world
        aExtent = [dLongitude_center -90, dLongitude_center + 90, -90,  90]
        minx, maxx, miny, maxy = aExtent

    if dLinewidth_in is not None:
        dLinewidth = dLinewidth_in
    else:
        dLinewidth = 0.25

    if sTitle_in is not None:
        sTitle = sTitle_in
        iFlag_title = 1
    else:
        iFlag_title = 0
        sTitle = ""

    if sFont_in is not None:
        sFont = sFont_in
    else:
        sFont = "Times New Roman"

    if sFilename_output_in is None:
        plt.ion()

    plt.rcParams["font.family"] = "DeJavu Serif"
    plt.rcParams["font.serif"] = sFont
    plt.rcParams["mathtext.fontset"] = "dejavuserif"

    fig = plt.figure(dpi=iDPI)
    fig.set_figwidth(iSize_x)
    fig.set_figheight(iSize_y)


    pDatasets_in = nc.Dataset(sFilename_mpas_mesh_netcdf_in)
    netcdf_format = pDatasets_in.file_format
    # read new netcdf
    for sKey, aValue in pDatasets_in.variables.items():
        # we need to filter out unused grids based on mpas specs
        if sKey == "latCell":
            latCell0 = aValue
        else:
            pass
        if sKey == "lonCell":
            lonCell0 = aValue
        else:
            pass

        if sKey == "edgesOnCell":
            edgesOnCell0 = aValue
        else:
            pass

        if sKey == "cellsOnCell":
            cellsOnCell0 = aValue
        else:
            pass

        if sKey == "verticesOnCell":
            verticesOnCell0 = aValue
        else:
            pass

        if sKey == "verticesOnEdge":
            verticesOnEdge0 = aValue
        else:
            pass

        if sKey == "indexToCellID":
            indexToCellID0 = aValue
        else:
            pass

        if sKey == "lonVertex":
            lonVertex0 = aValue
        else:
            pass

        if sKey == "latVertex":
            latVertex0 = aValue
        else:
            pass

    aLatitudeVertex = latVertex0[:] / math.pi * 180
    aLongitudeVertex = lonVertex0[:] / math.pi * 180
    # convert unit
    aLatitudeCell = latCell0[:] / math.pi * 180
    aLongitudeCell = lonCell0[:] / math.pi * 180
    aCellsOnCell = cellsOnCell0[:]
    # aCellOnEdge = cellsOnEdge0[:]
    aEdgesOnCell = edgesOnCell0[:]
    aVertexOnCell = verticesOnCell0[:]
    aVertexOnEdge0 = verticesOnEdge0[:]
    aIndexToCellID = indexToCellID0[:]

    aLongitudeCell_180 = convert_360_to_180(aLongitudeCell)
    aLongitudeVertex_180 = convert_360_to_180(aLongitudeVertex)

    if pProjection_map_in is not None:
        pProjection_map = pProjection_map_in
    else:
        dLon_mean = dLongitude_center
        dLat_mean = dLatitude_center
        pProjection_map = ccrs.Orthographic(
            central_longitude=dLon_mean, central_latitude=dLat_mean
        )
        print(dLon_mean, dLat_mean)

    if pProjection_data_in is not None:
        pProjection_data = pProjection_data_in
    else:
        pProjection_data = pSRS_wgs84

    pProjection_map._threshold /= 1.0e6

    ax = fig.add_axes([0.08, 0.1, 0.62, 0.7], projection=pProjection_map)
    plot_width_inch = fig.get_size_inches()[0] * fig.dpi
    char_width_inch = 0.1 * fig.dpi
    cwidth = int(plot_width_inch / char_width_inch)

    try:
        dAlpha = 1.0
        # only one of the base map can be used
        if iBasemap_zoom_level_in is not None:
            iBasemap_zoom_level = iBasemap_zoom_level_in
        else:
            image_size = [1000, 1000]
            scale_denominator = RasterTileServer.calculate_scale_denominator(aExtent, image_size)
            pSrc = osr.SpatialReference()
            pSrc.ImportFromEPSG(3857)  # mercator
            pProjection = pSrc.ExportToWkt()
            iBasemap_zoom_level = RasterTileServer.calculate_zoom_level(
                scale_denominator, pProjection, dpi=int(iDPI)
            )
            print('Basemap zoom level: ',iBasemap_zoom_level)
            pass
        if aBasemap_provider_in is not None:
            nTile_provider = len(aBasemap_provider_in)
            if nTile_provider > 5:
                print("Too many tile providers, only the first 5 will be used")
                nTile_provider = 5
            dAlpha = 1.0
            aLicense_info_list = []  # Collect all license info
            #for sBasemap_provider in aBasemap_provider_in:
            for i in range(nTile_provider):
                sBasemap_provider = aBasemap_provider_in[i]
                # Check if this is a Cartopy built-in provider (e.g., 'OSM')
                if sBasemap_provider == 'OSM':
                    # Use Cartopy's built-in OSM directly
                    from cartopy.io.img_tiles import OSM
                    pTiles = OSM()
                    ax.add_image(pTiles, iBasemap_zoom_level, alpha=dAlpha - i * 0.1)
                    # Add OSM attribution
                    aLicense_info_list.append("© OpenStreetMap contributors")
                else:
                    # Use RasterTileServer for custom providers
                    pTile_service = RasterTileServer(sBasemap_provider)
                    pTiles = pTile_service.get_cartopy_source_with_supersample(supersample = 2)
                    # For supersample, use higher zoom level and interpolation
                    ax.add_image(pTiles, iBasemap_zoom_level, alpha=dAlpha - i * 0.1, interpolation='bilinear')
                    aLicense_info_list.append(pTile_service.get_license_info())

            # Combine all license info and display once
            if aLicense_info_list:
                sLicense_info = " | ".join(aLicense_info_list)
                sLicense_info_wrapped = "\n".join(
                    textwrap.wrap(sLicense_info, width=cwidth)
                )
                ax.text(0.5,
                    0.05,
                    sLicense_info_wrapped,
                    transform=ax.transAxes,
                    ha="center",
                    va="center",
                    fontsize=6,
                    color="gray",
                    bbox=dict(
                        facecolor="white", edgecolor="black", boxstyle="round,pad=0.3"
                    ),
                )
            dAlpha = 0.5

    except URLError as e:
        print("No internet connection")
        dAlpha = 1.0

    nCell = len(aIndexToCellID)

    # Pre-filter cells based on extent to speed up plotting
    # Add buffer to extent to catch cells that might overlap the boundary
    extent_buffer = max((maxx - minx) * 0.05, (maxy - miny) * 0.05)
    minx_buffered = minx - extent_buffer
    maxx_buffered = maxx + extent_buffer
    miny_buffered = miny - extent_buffer
    maxy_buffered = maxy + extent_buffer

    # Find cells within the extent (using cell centers as initial filter)
    aCell_mask = np.logical_and(
        np.logical_and(aLongitudeCell_180 >= minx_buffered, aLongitudeCell_180 <= maxx_buffered),
        np.logical_and(aLatitudeCell >= miny_buffered, aLatitudeCell <= maxy_buffered)
    )
    aCell_indices = np.where(aCell_mask)[0]

    print(f"Processing {len(aCell_indices)} cells out of {nCell} total cells within extent")

    aPolygon = list()
    aEdgecolor = list()
    for i in aCell_indices:
        #dLon_center = float(aLongitudeCell_180[i])
        #dLat_center = float(aLatitudeCell[i])
        #lCellID = int(aIndexToCellID[i])
        aVertexOnCellIndex = np.array(aVertexOnCell[i,:])
        dummy0 = np.where(aVertexOnCellIndex > 0)
        aVertexIndex = aVertexOnCellIndex[dummy0] - 1
        aLonVertex = aLongitudeVertex_180[aVertexIndex]
        aLatVertex = aLatitudeVertex[aVertexIndex]
        nVertex = len(aLonVertex)
        if nVertex < 3:
            print("Vertex number is: ", nVertex, i)
        else:
            aCoords_gcs = list(zip(aLonVertex, aLatVertex))
            #do i need to close the polygon?  yes
            aCoords_gcs.append(aCoords_gcs[0])
            aPolygon.append(aCoords_gcs)
            aEdgecolor.append("none")

    if sColor_in is not None:
        sColor = sColor_in
    else:
        sColor = "black"
    if iFlag_fill == True:
        aPatch = [
            mpolygon(poly, closed=True, fill=iFlag_fill) for poly in aPolygon
        ]
        pPC = PatchCollection(
            aPatch,
            alpha=dAlpha,
            edgecolor="none",
            facecolor=sColor,
            transform=pProjection_data,
        )
    else:
        aPatch = [
            mpolygon(poly, closed=True, fill=iFlag_fill) for poly in aPolygon
        ]
        pPC = PatchCollection(
            aPatch,
            alpha=dAlpha,
            edgecolor=sColor,
            facecolor="none",
            linewidths=dLinewidth,
            transform=pProjection_data,
        )

    ax.add_collection(pPC)
    ax.set_extent(aExtent, crs=pSRS_wgs84)
    gl = ax.gridlines(
        crs=ccrs.PlateCarree(),
        draw_labels=True,
        linewidth=1,
        color="gray",
        alpha=0.5,
        linestyle="--",
        xlocs=np.arange(minx, maxx + (maxx - minx) / 9, (maxx - minx) / 8),
        ylocs=np.arange(miny, maxy + (maxy - miny) / 9, (maxy - miny) / 8),
    )
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlocator = mpl.ticker.MaxNLocator(4)
    gl.ylocator = mpl.ticker.MaxNLocator(4)
    gl.xlabel_style = {"size": 10, "color": "k", "rotation": 0, "ha": "right"}
    gl.ylabel_style = {"size": 10, "color": "k", "rotation": 90, "weight": "normal"}

    if iFlag_zebra == 1:
        ax.set_xticks(np.arange(minx, maxx + (maxx - minx) / 11, (maxx - minx) / 10))
        dummy = (maxy - miny) / 10
        ax.set_yticks(np.arange(miny, maxy + (maxy - miny) / 11, (maxy - miny) / 10))
        ax.set_axis_off()

    sTitle = "\n".join(textwrap.wrap(sTitle, width=cwidth))
    if iFlag_title is None:
        ax.set_title(sTitle)
    else:
        if iFlag_title == 1:
            ax.set_title(sTitle)
        else:
            pass

    if aLegend_in is not None:
        nlegend = len(aLegend_in)
        dLocation0 = 0.96
        for i in range(nlegend):
            sText = aLegend_in[i]
            dLocation = dLocation0 - i * 0.06
            ax.text(
                0.05,
                dLocation,
                sText,
                verticalalignment="top",
                horizontalalignment="left",
                transform=ax.transAxes,
                color="black",
                fontsize=iFont_size + 2,
            )

    if iFlag_zebra == 1:
        # ax.set_axis_off()
        ax.set_extent(aExtent, crs=pSRS_wgs84)
        ax.zebra_frame(crs=pSRS_wgs84, iFlag_outer_frame_in=1)

    ax.set_extent(aExtent, crs=pSRS_wgs84)
    if sFilename_output_in is None:
        plt.show(block=True)
        print("Finished plotting interactive map")
    else:
        if os.path.exists(sFilename_output_in):
            os.remove(sFilename_output_in)

        sDirname = os.path.dirname(sFilename_output_in)
        sFilename = os.path.basename(sFilename_output_in)
        sFilename_out = os.path.join(sDirname, sFilename)
        sExtension = os.path.splitext(sFilename)[1]
        if sExtension == ".png" or sExtension == ".jpg" or sExtension == ".jpeg":
            plt.savefig(sFilename_out, bbox_inches="tight")
        else:
            if sExtension == ".pdf":
                plt.savefig(sFilename_out, bbox_inches="tight")
            else:
                plt.savefig(sFilename_out, bbox_inches="tight", format="ps")

        plt.close("all")
        plt.clf()
        print("The plot is saved to: ", sFilename_out)