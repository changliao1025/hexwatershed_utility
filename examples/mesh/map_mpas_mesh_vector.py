from pyearthviz.map.vector.map_multiple_vector_files import map_multiple_vector_files

aFiletype_in = [3, 1]
aFilename_in = []

sFilename_mesh = "/Users/changliao/scratch/earthmesh/huanghe/mpas20261101002/mpas.geojson"
sFilename_river_raw = "/Users/changliao/scratch/earthmesh/huanghe/river_network/3.00E+03_6.40E+09/cropped_river_network_simplified_3.00E+03_6.40E+09.geojson"
aFilename_in.append(sFilename_mesh)
aFilename_in.append(sFilename_river_raw)

sFilename_out = "/Users/changliao/scratch/earthmesh/huanghe/mpas20261101002/mpas_river_map.png"

map_multiple_vector_files(
    aFiletype_in,
    aFilename_in,
    iFlag_title_in=None,
    iFlag_coastline_in=None,
    iFlag_zebra_in=1,
    iFlag_filter_in=None,
    iFlag_arrow_in=None,
    aFlag_thickness_in=[1, 1],
    aFlag_color_in=[1, 1],
    aFlag_colorbar_in=None,
    aFlag_discrete_in=None,
    aFlag_fill_in=None,
    aVariable_in=None,
    sFilename_output_in=sFilename_out,
    iFlag_scientific_notation_colorbar_in=None,
    aBasemap_provider_in=["Esri.NatGeo"],
    iFont_size_in=None,
    iBasemap_zoom_level_in=8,
    sColormap_in=None,
    sTitle_in=None,
    iDPI_in=600,
    iSize_x_in=None,
    iSize_y_in=None,
    aThickness_in=[0.3, 1],
    aColor_in=["black", "blue"],
    aMissing_value_in=None,
    aData_max_in=None,
    aData_min_in=None,
    sExtend_in=None,
    sFont_in=None,
    sUnit_in=None,
    aLegend_in=None,
    aExtent_in=None,
    pProjection_map_in=None,
    pProjection_data_in=None,
)
