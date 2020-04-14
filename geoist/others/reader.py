# Name        : reader.py
# Created on  : 2018/09/11 17:00
# Author      : Steve Chen<chenshi@cea-igp.ac.cn>
# Affiliation : Institute of Geophysics, CEA.
# Version     : 0.1.0
# License     : Distributed under the MIT License. See LICENSE.txt for more info.
# Github      : https://igp-gravity.github.io/
# Description : This file is part of GEOIST, which is free software for geophysicist
#               who need to process and analysis data and study inversion problem 
#               and distributed in the hope that it will be useful. Anyone can 
#               download and know news from above Github website address.     
#          
# Copyright (C) 2018-2020 GEOIST Development Team. All Rights Reserved.

# -*- coding: utf-8 -*-
"""
 Description : Application for grid processing.
"""

# third party imports
import numpy as np
import rasterio
import h5py

# local imports
from .grid2d import Grid2D
from .geodict import GeoDict


def _is_hdf(filename):
    """Detect whether file is a HDF format.

    Args:
        filename (str): Filename in question.
    Returns:
        bool: True if HDF, False if not.
    """
    is_hdf = False
    with open(filename, 'rb') as f:
        f.seek(1)
        bytes = header = f.read(3)
        try:
            header = bytes.decode()
            if header == 'HDF':
                is_hdf = True
        except UnicodeDecodeError:
            pass
        except Exception as e:
            raise e
    return is_hdf


def _get_geodict_from_src(src):
    """Return a geodict from a dataset object.

    Args:
        src (DatasetReader): Open rasterio DatasetReader object.
    Returns:
        GeoDict: GeoDict describing the DatasetReader object.
    """
    affine = src.transform
    geodict = {}
    geodict['dx'] = affine.a
    geodict['dy'] = -1 * affine.e
    geodict['xmin'] = affine.xoff + geodict['dx'] / 2.0
    geodict['ymax'] = affine.yoff - geodict['dy'] / 2.0
    geodict['ny'] = src.height
    geodict['nx'] = src.width
    geodict['xmax'] = geodict['xmin'] + (geodict['nx'] - 1) * geodict['dx']
    geodict['ymin'] = geodict['ymax'] - (geodict['ny'] - 1) * geodict['dy']

    gd = GeoDict(geodict)
    return gd


def _get_geodict_from_window(affine, window, data):
    """Return a geodict from a rasterio Window object.

    Args:
        affine (Affine): Object describing relationship between pixels
                         and coordinates.
        window (Window): Object describing the subset of the data that
                         was read from the file.
        data (ndarray): Data which was read from the file.
    Returns:
        GeoDict: GeoDict describing the Window of data read in.
    """
    xmin, ymax = affine * (window.col_off, window.row_off)
    geodict = {}
    geodict['dx'] = affine.a
    geodict['dy'] = -1 * affine.e
    geodict['xmin'] = xmin + geodict['dx'] / 2.0
    geodict['ymax'] = ymax - geodict['dy'] / 2.0
    nrows, ncols = data.shape
    geodict['ny'] = nrows
    geodict['nx'] = ncols
    geodict['xmax'] = geodict['xmin'] + (geodict['nx'] - 1) * geodict['dx']
    geodict['ymin'] = geodict['ymax'] - (geodict['ny'] - 1) * geodict['dy']

    gd = GeoDict(geodict)
    return gd


def _geodict_to_window(geodict, src, pad=False):
    """Convert a GeoDict to a rasterio Window object.

    Args:
        geodict (GeoDict): GeoDict describing the subset we wish to read.
        src (DatasetReader): Open rasterio DatasetReader object.
        pad (bool): Should the row/column indices be padded to include more
                    data?
    Returns:
        Window: Object describing the subset of the data that
                we wish to read from the file.
    """
    affine = src.transform
    dx = geodict.dx
    dy = geodict.dy

    # if padding requested, add one set of pixels all the way around
    xmin = geodict.xmin
    xmax = geodict.xmax
    ymin = geodict.ymin
    ymax = geodict.ymax
    if pad:
        xmin -= dx
        xmax += dx
        ymin -= dy
        ymax += dy

    # convert from pixel registered to gridline registered
    west = xmin - dx / 2.0
    south = (ymin + dy / 2.0) - dy
    east = (xmax - dx / 2.0) + dx
    north = ymax + dy / 2.0

    # handle meridian crossing by setting maximum to right edge of grid
    if east < west:
        east, _ = affine * (src.width, 0)

    # clip bounds to edges of grid file
    src_xmin, src_ymax = affine * (0, 0)
    src_xmax, src_ymin = affine * (src.width, src.height)

    if west < src_xmin:
        west = src_xmin
    if east > src_xmax:
        east = src_xmax
    if south < src_ymin:
        south = src_ymin
    if north > src_ymax:
        north = src_ymax

    window = src.window(west, south, east, north)
    row_range, col_range = window.toranges()

    # due to floating point rounding errors,
    # we may not get the integer row/col offsets we
    # actually want, so we're rounding here because
    # they're certain to be *very* close to the
    # integer values we want.
    row_start = int(np.round(row_range[0]))
    row_end = int(np.round(row_range[1]))
    col_start = int(np.round(col_range[0]))
    col_end = int(np.round(col_range[1]))
    bottom = max(0, row_start)
    top = min(src.height, row_end)
    row_range = (bottom, top)

    left = max(0, col_start)
    right = min(src.width, col_end)
    col_range = (left, right)

    window = rasterio.windows.Window.from_slices(row_range, col_range)
    return window


def _read_pixels(src, window):
    """Read pixels from a rasterio supported file format.

    NB: At the time of this writing, rasterio reading of
    compressed HDF files can be *very* slow, so for this format
    we are closing the input DatasetReader object, opening the
    file with h5py, reading the data, closing the file, and re-creating
    the DatasetReader object. This is why src is returned from the
    function.

    Args:
        src (DatasetReader): Open rasterio DatasetReader object.
        window (Window): Object describing the subset of the data that
                         we wish to read from the file.
    Returns:
        tuple: Numpy array containing data read, and a DatasetReader
               object representing the open file (see NB above).

    """
    fname = src.files[0]
    is_hdf = _is_hdf(fname)
    if src.driver != 'netCDF' or not is_hdf:
        data = np.squeeze(src.read(window=window), axis=0)
    else:
        src.close()
        f = h5py.File(fname, 'r')
        if window is None:
            data = np.flipud(f['z'])
        else:
            cstart = int(window.col_off)
            cend = cstart + int(window.width)
            rend = src.height - int(window.row_off)
            rstart = rend - int(window.height)
            # t1 = time.time()
            data = np.flipud(f['z'][rstart:rend, cstart:cend])
            # t2 = time.time()
            # print('h5py read: %.3f seconds.' % (t2-t1))
        f.close()
        src = rasterio.open(fname)

    # some kinds of files have NaN values encoded as special values like
    # -9999. I would have thought that rasterio w
    return (data, src)


def _read_data(src, samplegeodict, resample, method):
    """Read data from an open file, given subsetting/sampling information.

    This method will handle reading across the 180 meridian in the case
    of a global file, and a samplegeodict that spans that meridian.

    Args:
        src (DatasetReader): Open rasterio DatasetReader object.
        samplegeodict (GeoDict): GeoDict describing the subset we wish to read.
        resample (bool): True if resampling should be performed.
        method (str): One of ('nearest','linear').
    Returns:
        Grid2D: Object containing data and geospatial information.
    """
    affine = src.transform
    dx = affine.a
    dy = -1 * affine.e
    is_edge = samplegeodict.xmax < samplegeodict.xmin
    # read an extra row/column in every direction when resampling
    pad = False
    if resample:
        pad = True
    if not is_edge:
        window = _geodict_to_window(samplegeodict,
                                    src,
                                    pad=pad)
        data, src = _read_pixels(src, window)
        gd = _get_geodict_from_window(affine, window, data)
        gd.nodata = src.nodata
        grid = Grid2D(data, gd)
        return grid
    else:
        # split the windowing into two pieces - xmin to right edge
        # of global grid and left edge of global grid to xmax
        lxmin = samplegeodict.xmin
        # get right edge of grid, convert to pixel registration
        lxmax, _ = affine * (src.width, 0)
        lxmax += dx / 2.0
        ymin = samplegeodict.ymin
        ymax = samplegeodict.ymax
        sample_left = GeoDict.createDictFromBox(lxmin, lxmax,
                                                ymin, ymax, dx, dy)
        # get the left edge of the grid
        rxmin, _ = affine * (0, 0)
        # convert to pixel registered
        rxmin += dx / 2.0
        rxmax = samplegeodict.xmax
        sample_right = GeoDict.createDictFromBox(rxmin, rxmax,
                                                 ymin, ymax, dx, dy)
        left_window = _geodict_to_window(sample_left, src, pad=pad)

        left_block, src = _read_pixels(src, left_window)

        right_window = _geodict_to_window(sample_right, src, pad=pad)
        right_block, src = _read_pixels(src, right_window)

        left_gd = _get_geodict_from_window(affine,
                                           left_window,
                                           left_block)

        right_gd = _get_geodict_from_window(affine,
                                            right_window,
                                            right_block)
        data = np.concatenate((left_block, right_block), axis=1)
        nrows, ncols = data.shape
        geodict = {'xmin': left_gd.xmin,
                   'xmax': right_gd.xmax,
                   'ymin': ymin,
                   'ymax': ymax,
                   'dx': dx,
                   'dy': dy,
                   'nx': ncols,
                   'ny': nrows}
        gd = GeoDict(geodict)
        grid = Grid2D(data, gd)
        return grid


def get_file_geodict(filename):
    """Get the GeoDict describing the entire file.

    Args:
        filename (str): rasterio supported file format.
    Returns:
        GeoDict: GeoDict: GeoDict describing the entire file.
    """
    src = rasterio.open(filename)
    gd = _get_geodict_from_src(src)
    return gd


def read(filename, samplegeodict=None, resample=False,
         method='linear', doPadding=False, padValue=np.nan,
         apply_nan=True, force_cast=True):
    """Read part or all of a rasterio file, resampling and padding as necessary.

    If samplegeodict is not provided, then the entire file will be read.

    If samplegeodict and resample are provided, then the smallest subset of
    data containing the samplegeodict plus a 1 pixel wide ring of data
    will be read in, then that data will be resampled to the samplegeodict
    bounds/resolution.

    If doPadding is set to True, then pixels on the edge of the source grid
    will be padded 1 pixel deep with input padValue.

    In addition to pad pixels, extra data may be read around the edges of
    the desired area to ensure raw data is read on integer row/column offsets.

    Args:
        filename (str): rasterio supported file format.
        samplegeodict (GeoDict): GeoDict describing the subset we wish to read.
        resample (bool): True if resampling should be performed.
        method (str): One of ('nearest','linear').
        doPadding (bool): Whether to add ring of padValue pixels after reading
                          from file.
        padValue (float): Value to insert in ring of padding pixels.
        apply_nan (bool): Convert nodata values to NaNs, upcasting to float if necessary.
        force_cast (bool): If data values exceed 
    """
    # use rasterio to read all formats
    src = rasterio.open(filename)

    # first establish if we are subsetting the data.
    # if not, read the whole file and return.
    if samplegeodict is None:
        data, src = _read_pixels(src, None)
        gd = _get_geodict_from_src(src)
        grid = Grid2D(data, gd)
        src.close()
        return grid

    # if non-nearest resampling, this grid may have a ring of padding pixels
    # around the outside.
    grid = _read_data(src, samplegeodict, resample, method)

    # make sure this raw grid is big enough to support resampling
    # if padding is turned off
    rdict = grid.getGeoDict()
    if resample:
        if method == 'nearest':
            c1 = rdict.xmin <= samplegeodict.xmin
            c2 = rdict.xmax >= samplegeodict.xmax
            c3 = rdict.ymin <= samplegeodict.ymin
            c4 = rdict.ymax >= samplegeodict.ymax
        else:
            dx = rdict.dx
            dy = rdict.dy
            c1 = rdict.xmin > samplegeodict.xmin - dx
            c2 = rdict.xmax < samplegeodict.xmax + dx
            c3 = rdict.ymin > samplegeodict.ymin - dy
            c4 = rdict.ymax < samplegeodict.ymax + dy
        if (c1 or c2 or c3 or c4) and not doPadding:
            msg = ('Without padding you cannot ask to resample from '
                   'edge of grid.')
            raise IndexError(msg)

    grid._geodict.nodata = src.nodata
    if apply_nan:
        grid.applyNaN(force=force_cast)

    if doPadding:
        filedict = get_file_geodict(filename)
        # use the padDict method of Grid2D to create our padded grid
        # Pad one row/col on all sides.
        pd = grid.getPadding(filedict, samplegeodict, doPadding=True)
        data, gd = Grid2D.padGrid(grid._data, grid._geodict, pd)
        if len(data[np.isinf(data)]):
            data[np.isinf(data)] = padValue
        grid = Grid2D(data, gd)

    if resample:
        grid = grid.interpolateToGrid(samplegeodict, method=method)

    src.close()
    return grid
