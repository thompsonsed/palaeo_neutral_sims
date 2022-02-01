# File containing the functions relating to the data preparation of paleo maps and sampling points
# Here we want to
# - convert the shape files to rasters
# - create the sample mask
# - record the sites of the sample points (with each site referred to as a 'fragment' within the codebase)
# - record the times of each sample point
import random

from osgeo import gdal, ogr, osr
import os
import math
import pandas as pd
import numpy as np
from pycoalescence import Map


def extract_from_shapefile(src_file, dest_file, field, field_value):
    """
    Extracts elements from the source file if the field matches the field value.
    Writes the output to a new shape file.

    :note This will delete any existing destination file

    :param src_file: the source shape file to open
    :param dest_file: the path to the destination shape file
    :param field: the field to check against
    :param field_value: the value in the field to extract values from

    :rtype None
    """
    if not src_file.endswith(".shp"):
        raise ValueError("Input file {} is not a shape file.".format(src_file))
    if not os.path.exists(src_file):
        raise IOError("Input file {} does not exist.".format(src_file))
    if os.path.exists(dest_file):
        os.remove(dest_file)
    # Get the input Layer
    src_driver = ogr.GetDriverByName("ESRI Shapefile")
    src_ds = src_driver.Open(src_file, 0)
    if src_ds is None:
        raise IOError("Could not open {}".format(src_file))
    src_layer = src_ds.GetLayer()
    src_layer.SetAttributeFilter("{} = '{}'".format(field, field_value))
    # Create the destination layers
    dest_driver = ogr.GetDriverByName("ESRI Shapefile")
    # Create the destination shapefile
    dest_ds = dest_driver.CreateDataSource(dest_file)
    dest_layer_name = os.path.splitext(os.path.split(dest_file)[1])[0]
    dest_layer = dest_ds.CreateLayer(dest_layer_name, geom_type=ogr.wkbMultiPolygon)
    # Add input Layer Fields to the destination layer
    src_layer_defn = src_layer.GetLayerDefn()
    for i in range(0, src_layer_defn.GetFieldCount()):
        field_defn = src_layer_defn.GetFieldDefn(i)
        dest_layer.CreateField(field_defn)
    # Get the output Layer's Feature Definition
    dest_layer_defn = dest_layer.GetLayerDefn()
    # Add features to the ouput Layer
    for src_feature in src_layer:
        # Create output Feature
        dest_feature = ogr.Feature(dest_layer_defn)
        # Add field values from input Layer
        for i in range(0, dest_layer_defn.GetFieldCount()):
            dest_feature.SetField(dest_layer_defn.GetFieldDefn(i).GetNameRef(), src_feature.GetField(i))
        # Set geometry as centroid
        geom = src_feature.GetGeometryRef()
        if geom is None:
            continue
        dest_feature.SetGeometry(geom.Clone())
        # Add new feature to output Layer
        dest_layer.CreateFeature(dest_feature)
        dest_feature = None
    # Save and close DataSources
    src_ds = None
    dest_ds = None


def find_distance_lat_long(start_lat, start_long, end_lat, end_long, res):
    """
    Finds the distance in number of cells at the resolution provided.

    :param start_lat:
    :param start_long:
    :param end_lat:
    :param end_long:
    :param res:
    :return:
    """
    dist_lat = abs(start_lat - end_lat) / res
    dist_long = abs(start_long - end_long) / res
    return [dist_lat, dist_long]


def rasterize(vector_file, raster_file, pixel_size=1, EPSG=4326, field=None, default_val=0, dataType=gdal.GDT_Float32):
    """Converts the vector file to a raster file.

    :note Note this will delete any existing raster file.

    :param vector_file: path to the shapefile for conversion
    :param raster_file: path to the output raster file to be created
    :param pixel_size: optionally provide the pixel size
    :param field: optionally, the name of the field to set as the raster value
    :param default_val: the no data value

    :rtype None
    """
    # Only continue if the input file is a valid shape file and it exists.
    if not vector_file.endswith(".shp"):
        raise ValueError("Input file {} is not a shape file.".format(vector_file))
    if not os.path.exists(vector_file):
        raise IOError("Input file {} does not exist.".format(vector_file))
    # Delete the raster file, if it currently exists.
    if os.path.exists(raster_file):
        os.remove(raster_file)
    # Open the vector file
    orig_data_src = ogr.Open(vector_file)
    # Make a copy of the layer's data source because we'll need to
    # modify its attributes table
    source_ds = ogr.GetDriverByName("Memory").CopyDataSource(orig_data_src, "")
    source_layer = source_ds.GetLayer(0)
    source_srs = source_layer.GetSpatialRef()
    x_min, x_max, y_min, y_max = source_layer.GetExtent()
    x_min -= pixel_size * 5.5
    y_max += pixel_size * 4.5
    x_res = math.ceil((x_max - x_min) / pixel_size) + 10
    y_res = math.ceil((y_max - y_min) / pixel_size) + 10
    target_ds = gdal.GetDriverByName("GTiff").Create(raster_file, x_res, y_res, 1, dataType)
    target_ds.SetGeoTransform((x_min, pixel_size, 0, y_max, 0, -pixel_size))
    if source_srs:
        # Make the target raster have the same projection as the source
        target_ds.SetProjection(source_srs.ExportToWkt())
    else:
        # Source has no projection (needs GDAL >= 1.7.0 to work)
        target_ds.SetProjection('LOCAL_CS["arbitrary"]')
    # Rasterize
    if field is None:
        # print('None')
        err = gdal.RasterizeLayer(target_ds, [1], source_layer, None, None, [1], ["ALL_TOUCHED=FALSE"])
    else:
        err = gdal.RasterizeLayer(
            target_ds,
            [1],
            source_layer,
            options=["ALL_TOUCHED=FALSE", "ATTRIBUTE=" + field, "N0_DATA=" + str(default_val)],
        )
    if err != 0:
        raise IOError("Gdal error while rasterising layer. Error code: {}".format(err))


def get_geo_info(src_file):
    """
    Gets the geo info data for a raster file.
    :return: dictionary containing the no_data_value, x_res, y_res, geo_transform, projection and data_type
    :rtype dict
    """
    if not os.path.exists(src_file):
        raise IOError("Source file {} does not exist.".format(src_file))
    src_ds = gdal.Open(src_file, gdal.GA_ReadOnly)
    src_band = src_ds.GetRasterBand(1)
    NDV = src_band.GetNoDataValue()
    x_res = src_ds.RasterXSize
    y_res = src_ds.RasterYSize
    geo_transform = src_ds.GetGeoTransform()
    projection = osr.SpatialReference()
    projection.ImportFromWkt(src_ds.GetProjectionRef())
    data_type = gdal.GetDataTypeName(src_band.DataType)
    src_ds = None
    return {
        "no_data_value": NDV,
        "x_res": x_res,
        "y_res": y_res,
        "geo_transform": geo_transform,
        "data_type": data_type,
    }


def create_shapefile_from_points(data_frame: pd.DataFrame, dest_file):
    """
    Creates a shapefile from a pandas dataframe with a lat and long column
    Creates one point per line in the data frame.

    :param data_frame: the pandas dataframe, containing a lat and long column
    :param dest_file: the location of the output shape file
    """
    # Create the destination layers
    dest_driver = ogr.GetDriverByName("ESRI Shapefile")
    # Create the destination shapefile
    dest_ds = dest_driver.CreateDataSource(dest_file)
    dest_layer_name = os.path.splitext(os.path.split(dest_file)[1])[0]
    dest_layer = dest_ds.CreateLayer(dest_layer_name, geom_type=ogr.wkbPoint)
    # Add input Layer Fields to the destination layer
    fields_list = data_frame.columns.values
    for field_defn in fields_list:
        field_ref = ogr.FieldDefn(field_defn, ogr.OFTReal)
        if "Unnamed" in field_defn:
            continue
        dest_layer.CreateField(field_ref)
    # Get the output Layer's Feature Definition
    dest_layer_defn = dest_layer.GetLayerDefn()
    # Add features to the ouput Layer
    for index, row in data_frame.iterrows():
        # Create output Feature
        dest_feature = ogr.Feature(dest_layer_defn)
        # Add field values from input Layer
        for field in fields_list:
            if "lat" != field and "long" != field:
                dest_feature.SetField(field, row[field])
        # Set geometry as centroid
        point = ogr.Geometry(ogr.wkbPoint)
        point.AddPoint(float(row["long"]), float(row["lat"]))
        dest_feature.SetGeometry(point)
        # Add new feature to output Layer
        dest_layer.CreateFeature(dest_feature)
        dest_feature = None
    # Save and close DataSources
    src_ds = None
    dest_ds = None


def mask_value(input_file, output_file, mask_val):
    """
    Generates a mask raster for the input raster matching the mask value in a new file, the output raster.

    :param input_file: the input file to select values from
    :param output_file: the output file location
    :param mask_val: the value to use as the mask value

    :rtype None
    """
    src_ds = gdal.Open(input_file, gdal.GA_ReadOnly)
    src_band = src_ds.GetRasterBand(1)
    raster_arr = np.array(src_band.ReadAsArray())
    # Create a new file for outputting t
    x_res = src_ds.RasterXSize
    y_res = src_ds.RasterYSize
    target_ds = gdal.GetDriverByName("GTiff").Create(output_file, x_res, y_res, 1, dataType)
    target_ds.SetGeoTransform((x_min, pixel_size, 0, y_max, 0, -pixel_size))
    mask = raster_arr == mask_val


def randomly_clear_landscape(input_file, output_file, proportion_cover):
    """
    Randomly clears the landscape so that habitat cover present in input_file is reduced to proportion_cover in
    output_file.
    :param input_file: the file to read the landscape from
    :param output_file: the output path to store the file to
    :param proportion_cover: the proportion of cover to reduce the habitat to
    """
    ds = gdal.Open(input_file)
    band = ds.GetRasterBand(1)
    arr = band.ReadAsArray()
    [cols, rows] = arr.shape
    # Modify the array based on randomly removing pixels
    total = np.sum(arr)
    desired_total = math.floor(total * proportion_cover)
    to_remove = total - desired_total
    # input array
    temp = arr.flatten()  # Flatten to 1D
    inds = np.random.choice(2, size=temp.size, p=[proportion_cover, 1 - proportion_cover]).astype(np.bool)
    # Create masked array of inds
    temp = np.ma.array(temp, mask=inds, fill_value=0).filled()
    temp = temp.reshape(arr.shape)
    driver = gdal.GetDriverByName("GTiff")
    geo_info = get_geo_info(input_file)
    outdata = driver.Create(output_file, rows, cols, 1, gdal.GDT_Int16)
    outdata.SetGeoTransform(ds.GetGeoTransform())  ##sets same geotransform as input
    outdata.SetProjection(ds.GetProjection())  ##sets same projection as input
    outdata.GetRasterBand(1).WriteArray(temp)
    outdata.GetRasterBand(1).SetNoDataValue(10000)  ##if you want these values transparent
    outdata.FlushCache()  ##saves to disk!!
    outdata = None
    band = None
    ds = None


def add_masks(input_file_a, input_file_b, output_file):
    """
    Adds together the two input files and combines them into the output file.
    Using the np.logical_or function, so if the mask is true in either location, the output will be true.
    :param input_file_a: the first file to add, the larger of the two if they're not the same size
    :param input_file_b: the second file to add, the smaller
    :param output_file: the path to the output file
    """
    ds_a = gdal.Open(input_file_a)
    band_a = ds_a.GetRasterBand(1)
    arr_a = band_a.ReadAsArray()
    [rows_a, cols_a] = arr_a.shape
    ds_b = gdal.Open(input_file_b)
    band_b = ds_b.GetRasterBand(1)
    arr_b = band_b.ReadAsArray()
    [rows_b, cols_b] = arr_b.shape
    # check dimensions add up
    if cols_a < cols_b or rows_a < rows_b:
        raise ValueError(
            "First input file must be larger than second. ({}, {}) smaller than ({},{})".format(
                cols_a, rows_a, cols_b, rows_b
            )
        )
    if 0 in [cols_a, cols_b, rows_a, rows_b]:
        raise ValueError(
            "One of the files has a dimension of 0: ({}, {}), ({}, {})".format(cols_a, rows_a, cols_b, rows_b)
        )
    ulx_a, xres_a, _, uly_a, _, yres_a = ds_a.GetGeoTransform()
    ulx_b, xres_b, _, uly_b, _, yres_b = ds_b.GetGeoTransform()
    offset_cols = round((ulx_b - ulx_a) / xres_a)
    offset_rows = round((uly_b - uly_a) / yres_b)
    # modify values in our array
    replacement = arr_b + arr_a[offset_rows : offset_rows + rows_b, offset_cols : offset_cols + cols_b]
    out = np.array(arr_a)
    out[offset_rows : offset_rows + rows_b, offset_cols : offset_cols + cols_b] = replacement
    driver = gdal.GetDriverByName("GTiff")
    outdata = driver.Create(output_file, cols_a, rows_a, 1, gdal.GDT_Int16)
    outdata.SetGeoTransform(ds_a.GetGeoTransform())  ##sets same geotransform as input
    outdata.SetProjection(ds_a.GetProjection())  ##sets same projection as input
    outdata.GetRasterBand(1).WriteArray(out)
    outdata.GetRasterBand(1).SetNoDataValue(10000)  ##if you want these values transparent
    outdata.FlushCache()
    # Clean up
    outdata = None
    ds_a = None
    ds_b = None
    band_a = None
    band_b = None


def create_cluster(input_map, output_map, coordinates, radius):
    """
    Creates clusters on the input map, setting values to 1 within radius of the given coordinates.

    :param str input_map: the input file path
    :param output_map: the output file path
    :param list coordinates: list of x, y coordinates to apply
    :param int radius: number of cell widths to include
    """
    if not os.path.exists(input_map):
        raise IOError("File does not exist at {}.".format(input_map))
    if os.path.exists(output_map):
        raise FileExistsError("File already exists at {}.".format(output_map))
    m_input = Map(input_map)
    m_input.create_copy(output_map)
    m_output = Map(output_map)
    dim_x, dim_y = [abs(x) for x in m_input.get_x_y()]
    m_output.data = np.zeros(shape=(dim_y, dim_x))
    m_output.band_number = 1
    m_output.write()
    m_output.data = None
    for x, y in coordinates:
        if x > dim_x or x < 0 or y > dim_y or y < 0:
            raise ValueError("Coordinates {}, {} are not in map of dimensions {}, {}.".format(x, y, dim_x, dim_y))
        x_min = x - radius
        x_max = x + radius
        if x_min < 0:
            x_min = 0
        if x_max > dim_x:
            x_max = dim_x
        x_size = x_max - x_min

        y_min = y - radius
        y_max = y + radius
        if y_min < 0:
            y_min = 0
        if y_max > dim_y:
            y_max = dim_y
        y_size = y_max - y_min
        tmp = np.zeros(shape=(y_size, x_size))
        extraction = m_input.get_subset(x_offset=x_min, y_offset=y_min, x_size=x_size, y_size=y_size)
        ys, xs = np.ogrid[-radius:radius, -radius:radius]
        index = xs ** 2 + ys ** 2 <= radius ** 2
        index = index[max(radius - y - y_min, 0) : y_max - y + radius, max(radius - x - x_min, 0) : x_max - x + radius]
        tmp[index] = 1
        tmp[extraction == 0] = 0
        extraction_output = m_output.get_subset(x_offset=x_min, y_offset=y_min, x_size=x_size, y_size=y_size)
        tmp = np.maximum(extraction_output, tmp)
        m_output.write_subset(tmp, x_off=x_min, y_off=y_min)
