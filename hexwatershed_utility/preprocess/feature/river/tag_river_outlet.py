import os, sys
from osgeo import gdal, ogr
import numpy as np
import logging
from pyearth.gis.gdal.read.raster.gdal_read_geotiff_file import gdal_read_geotiff_file
from pyearth.gis.gdal.write.raster.gdal_write_geotiff_file import gdal_write_geotiff_file
logger = logging.getLogger(__name__)
def tag_river_outlet(sFilename_flowline_hydrosheds_in,
                     sFilename_river_network_raster_in,
                     sFilename_river_network_raster_mouth_tagged,
                      iSpecial_value_outlet_in = 255):

    def _extract_last_vertex_xy(pGeometry_in):
        """Extract (x, y) of the last vertex from line-like geometry."""
        if pGeometry_in is None:
            return None
        sGeometry_name = pGeometry_in.GetGeometryName().upper()
        if sGeometry_name in ("LINESTRING", "LINEARRING"):
            nPoint_count = pGeometry_in.GetPointCount()
            if nPoint_count > 0:
                dX, dY, _ = pGeometry_in.GetPoint(nPoint_count - 1)
                return dX, dY
            return None
        if sGeometry_name == "POINT":
            dX = pGeometry_in.GetX()
            dY = pGeometry_in.GetY()
            return dX, dY
        if sGeometry_name in ("MULTILINESTRING", "MULTIPOINT", "GEOMETRYCOLLECTION"):
            nGeometry_count = pGeometry_in.GetGeometryCount()
            for iGeometry in range(nGeometry_count - 1, -1, -1):
                pSub_geometry = pGeometry_in.GetGeometryRef(iGeometry)
                pXY = _extract_last_vertex_xy(pSub_geometry)
                if pXY is not None:
                    return pXY
            return None
        return None


    #check input files
    if not os.path.isfile(sFilename_flowline_hydrosheds_in):
        raise FileNotFoundError(f"Input flowline file not found: {sFilename_flowline_hydrosheds_in}")
    if not os.path.isfile(sFilename_river_network_raster_in):
        raise FileNotFoundError(f"Input raster file not found: {sFilename_river_network_raster_in}")

    #remove output file if it exists
    if os.path.isfile(sFilename_river_network_raster_mouth_tagged):
        os.remove(sFilename_river_network_raster_mouth_tagged)

    #read input river network vector
    pDataset_vector = ogr.Open(sFilename_flowline_hydrosheds_in, 0)
    if pDataset_vector is None:
        raise RuntimeError(f"Failed to open vector file: {sFilename_flowline_hydrosheds_in}")
    pLayer_vector = pDataset_vector.GetLayer(0)
    if pLayer_vector is None:
        raise RuntimeError(f"Failed to get layer from vector file: {sFilename_flowline_hydrosheds_in}")

    #read input river network raster
    pDataset_dummy = gdal_read_geotiff_file(sFilename_river_network_raster_in)
    aData_river = pDataset_dummy['dataOut']
    dOrigin_x = pDataset_dummy['originX']
    dOrigin_y = pDataset_dummy['originY']
    dPixel_width = pDataset_dummy['pixelWidth']
    dPixel_height = pDataset_dummy['pixelHeight']
    nrow = pDataset_dummy['nrow']
    ncolumn = pDataset_dummy['ncolumn']
    dMissing_value_in = pDataset_dummy['missingValue']
    pProjection_in = pDataset_dummy['projection']
    dataType = pDataset_dummy['dataType']

    pLayer_vector.ResetReading()

    for pFeature in pLayer_vector:
        pOutlet_xy = None
        pGeometry = pFeature.GetGeometryRef()
        pOutlet_xy = _extract_last_vertex_xy(pGeometry)
        if pOutlet_xy is not None:
            iColumn = int(np.floor((pOutlet_xy[0] - dOrigin_x) / dPixel_width))
            iRow = int(np.floor((dOrigin_y - pOutlet_xy[1]) / abs(dPixel_height)))
            if 0 <= iRow < nrow and 0 <= iColumn < ncolumn:
                aData_river[iRow, iColumn] = iSpecial_value_outlet_in
                logger.info(f"Tagged outlet at ({iRow}, {iColumn}) with value {iSpecial_value_outlet_in}")

            else:
                logger.warning(
                    f"Outlet vertex ({pOutlet_xy[0]}, {pOutlet_xy[1]}) is out of raster bounds."
                )
        else:
            logger.warning(
                "iSpecial_value_outlet_in provided, but no line-like geometry with vertices was found"
            )

    #write output raster
    gdal_write_geotiff_file(
        sFilename_river_network_raster_mouth_tagged,
        aData_river,
        dPixel_width,
        dPixel_height,
        dOrigin_x,
        dOrigin_y,
        dMissing_value_in,
        pProjection_in,
        datatype=dataType,
    )

    print(f"Tagged river outlet raster written to: {sFilename_river_network_raster_mouth_tagged}")

