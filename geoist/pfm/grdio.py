# -*- coding: utf-8 -*-
"""
 Name        : grdio.py
 Created on  : 2018/11/24 08:57
 Author      : Steve Chen <chenshi@cea-igp.ac.cn>
 Affiliation : Institute of Geophysics, CEA.
 Version     : 0.1.0
 Copyright   : Copyright (C) 2018-2020 GEOIST Development Team. All Rights Reserved.
 License     : Distributed under the MIT License. See LICENSE.txt for more info.
 Github      : https://igp-gravity.github.io/
 Description : Application for processing grid data of potential dataset.
"""
import struct
import numpy as np
import scipy.interpolate as interp
from matplotlib import pyplot as plt
import warnings

def _check_area(area):
    """
    Check that the area argument is valid.
    For example, the west limit should not be greater than the east limit.
    """
    x1, x2, y1, y2 = area
    assert x1 <= x2, \
        "Invalid area dimensions {}, {}. x1 must be < x2.".format(x1, x2)
    assert y1 <= y2, \
        "Invalid area dimensions {}, {}. y1 must be < y2.".format(y1, y2)

def regular(area, shape, z=None):
    """
    Create a regular grid.

    The x directions is North-South and y East-West. Imagine the grid as a
    matrix with x varying in the lines and y in columns.

    Returned arrays will be flattened to 1D with ``numpy.ravel``.

    Parameters:

    * area
        ``(x1, x2, y1, y2)``: Borders of the grid
    * shape
        Shape of the regular grid, ie ``(nx, ny)``.
    * z
        Optional. z coordinate of the grid points. If given, will return an
        array with the value *z*.

    Returns:

    * ``[x, y]``
        Numpy arrays with the x and y coordinates of the grid points
    * ``[x, y, z]``
        If *z* given. Numpy arrays with the x, y, and z coordinates of the grid
        points

    Examples:

    >>> x, y = regular((0, 10, 0, 5), (5, 3))
    >>> print(x)
    [  0.    0.    0.    2.5   2.5   2.5   5.    5.    5.    7.5   7.5   7.5
      10.   10.   10. ]
    >>> print(x.reshape((5, 3)))
    [[  0.    0.    0. ]
     [  2.5   2.5   2.5]
     [  5.    5.    5. ]
     [  7.5   7.5   7.5]
     [ 10.   10.   10. ]]

    """
    nx, ny = shape
    x1, x2, y1, y2 = area
    _check_area(area)
    xs = np.linspace(x1, x2, nx)
    ys = np.linspace(y1, y2, ny)
    # Must pass ys, xs in this order because meshgrid uses the first argument
    # for the columns
    arrays = np.meshgrid(ys, xs)[::-1]
    if z is not None:
        arrays.append(z*np.ones(nx*ny, dtype=np.float))
    return [i.ravel() for i in arrays]

def spacing(area, shape):
    """
    Returns the spacing between grid nodes

    Parameters:

    * area
        ``(x1, x2, y1, y2)``: Borders of the grid
    * shape
        Shape of the regular grid, ie ``(nx, ny)``.

    Returns:

    * ``[dx, dy]``
        Spacing the y and x directions

    Examples:

    >>> print(spacing((0, 10, 0, 20), (11, 11)))
    [1.0, 2.0]
    >>> print(spacing((0, 10, 0, 20), (11, 21)))
    [1.0, 1.0]
    >>> print(spacing((0, 10, 0, 20), (5, 21)))
    [2.5, 1.0]
    >>> print(spacing((0, 10, 0, 20), (21, 21)))
    [0.5, 1.0]

    """
    x1, x2, y1, y2 = area
    nx, ny = shape
    dx = (x2 - x1)/(nx - 1)
    dy = (y2 - y1)/(ny - 1)
    return [dx, dy]

class grddata(object):
    """
    Grid Data Object
    Attributes
    ----------
    data : numpy masked array
        array to contain raster data
    xmin : float
        min value X coordinate of raster grid
    ymin : float
        min value Y coordinate of raster grid
    xdim : float
        x-dimension of grid cell
    ydim : float
        y-dimension of grid cell
    typeofdata : int
        number of datatype
    dataname : str
        data name or id
    rows : int
        number of rows for each raster grid/band
    cols : int
        number of columns for each raster grid/band
    nullvalue : float
        grid null or nodata value
    norm : dictionary
        normalized data
    gtr : tuple
        projection information
    wkt : str
        projection information
    units : str
        description of units to be used with color bars
    """
    def __init__(self):
        self.data = np.ma.array([])
        self.data0 = np.array([])
        self.xmin = 0.0  # min value of X coordinate
        self.ymin = 0.0  # min value of Y coordinate
        self.xdim = 1.0
        self.ydim = 1.0
        self.dmin = 0.0
        self.dmax = 0.0
        self.typeofdata = 1 # 1- grav or 2- mag
        self.dataname = '' #name of data
        self.rows = -1
        self.cols = -1
        self.nullvalue = 1e+20
        self.norm = {}
        self.gtr = (0.0, 1.0, 0.0, 0.0, -1.0)
        self.wkt = ''
        self.units = ''

    def fill_nulls(self, method='nearest'):
        """
            Fill in the NaNs or masked values on interpolated points using nearest
            neighbors.
            method='nearest' or 'linear' or 'cubic'
        """
        if np.ma.is_masked(self.data):
            nans = self.data.mask
        else:
            nans = np.isnan(self.data)

        nx,ny = nans.shape
        ns = nans.reshape(nx*ny)
        shape = (nx, ny)
        xmax = self.xmin + (self.cols-1)*self.xdim
        ymax = self.ymin + (self.rows-1)*self.ydim
        area = (self.xmin, xmax, self.ymin, ymax)
        x, y = regular(area, shape)
        dtmp = self.data.copy() #数组copy，不改变源数组
        dtmp1 = dtmp.reshape(nx*ny)
        ns1 = (ns == False)
        dtmp1[ns] = interp.griddata((x[ns1], y[ns1]), dtmp1[ns1], (x[ns], y[ns]),
                                    method).ravel()
        self.data0 = dtmp1.reshape(nx,ny)

    def grd2xyz(self, flag = True):
        """
        Return x,y,z 1-D array data from 2-D grid array.

        Parameters:
          flag  : True  -  Output Grid Grid
                False -  Output Bak Grid Grid
        Returns:
          x,y,z 1-D array data
        """
        nx,ny = self.data.shape
        xmax = self.xmin + (self.cols-1)*self.xdim
        ymax = self.ymin + (self.rows-1)*self.ydim

        shape = (nx, ny)
        area = (self.xmin, xmax, self.ymin, ymax)
        x, y = regular(area, shape)
        if flag:
          z = self.data.reshape(nx*ny)
        else:
          z = self.data0.reshape(nx*ny)
        return (x, y, z)


    def load_surfer(self, fname, *args, **kwargs):
        """
        Read data from a Surfer grid file.

        Parameters:

        * fname : str
            Name of the Surfer grid file
        * dtype : numpy dtype object or string
            The type of variable used for the data. Default is numpy.float64 for
            ascii data and is '=f' for binary data. Use numpy.float32 if the
            data are large and precision is not an issue.
        * header_format : header format (excluding the leading 'DSBB') following
            the convention of the struct module. Only used for binary data.

        Returns:

        """
        with open(fname,'rb') as f:
            tmp = f.read(4)
        if tmp == b'DSAA':
            self._load_surfer_ascii(fname,*args,**kwargs)
        elif tmp == b'DSBB':
            self._load_surfer_dsbb(fname,*args,**kwargs)
        else:
            raise ValueError("Unknown header info {}.".format(tmp)
                            +"Only DSAA or DSBB could be recognized.")

    def _load_surfer_dsbb(self,fname,dtype='=f',header_format='cccchhdddddd'):
        """
        Read data from a Surfer DSBB grid file.

        Parameters:

        * fname : str
            Name of the Surfer grid file
        * dtype : numpy dtype object or string
            The type of variable used for the data. Default is numpy.float64. Use
            numpy.float32 if the data are large and precision is not an issue.
        * header_format : header format following the convention of the
            struct module.

        Returns:

        """
        with open(fname,'rb') as f:
            # read header
            header_len = struct.calcsize(header_format)
            header = f.read(header_len)
            # read data
            data = b''
            for x in f:
                data += x

        # unpack header
        s = struct.Struct(header_format)
        (tmp,tmp,tmp,tmp,self.cols,self.rows,self.xmin,self.xmax,
         self.ymin,self.ymax,self.dmin,self.dmax) = s.unpack(header)
        if self.cols<=0 and self.rows<=0:
            raise ValueError("Array shape can't be infered.")

        # convert data to numpy array
        self.data = np.frombuffer(data,dtype=dtype).reshape(self.cols,self.rows)
        self.data = np.ma.MaskedArray(self.data)
        self.cols,self.rows = self.data.shape
        if self.data.min()+1<self.dmin or self.data.max()-1>self.dmax:
            warnings.warn("(min(z),max(z)) in the data is incompatible "
                          +"with (zmin,zmax) in the header. "
                          +"Please check whether the 'dtype' argument is "
                          +"correct.(default is '=f')")
        self.xdim = (self.xmax-self.xmin)/(self.rows-1)
        self.ydim = (self.ymax-self.ymin)/(self.cols-1)


    def _load_surfer_ascii(self, fname, dtype='float64'):
        """
        Read data from a Surfer ASCII grid file.

        Parameters:

        * fname : str
            Name of the Surfer grid file
        * dtype : numpy dtype object or string
            The type of variable used for the data. Default is numpy.float64. Use
            numpy.float32 if the data are large and precision is not an issue.

        Returns:

        """
        # Surfer ASCII grid structure
        # DSAA            Surfer ASCII GRD ID
        # nCols nRows     number of columns and rows
        # xMin xMax       X min max
        # yMin yMax       Y min max
        # zMin zMax       Z min max
        # z11 z21 z31 ... List of Z values
        with open(fname) as input_file:
            # DSAA is a Surfer ASCII GRD ID (discard it for now)
            input_file.readline()
            # Read the number of columns (ny) and rows (nx)
            ny, nx = [int(s) for s in input_file.readline().split()]
            #shape = (nx, ny)
            # Our x points North, so the first thing we read is y, not x.
            ymin, ymax = [float(s) for s in input_file.readline().split()]
            xmin, xmax = [float(s) for s in input_file.readline().split()]
            #area = (xmin, xmax, ymin, ymax)
            dmin, dmax = [float(s) for s in input_file.readline().split()]
            field = np.fromiter((float(s)
                                 for line in input_file
                                 for s in line.split()),
                                dtype=dtype)
            nans = field >= 1.70141e+38
            if np.any(nans):
                field = np.ma.masked_where(nans, field)
            #err_msg = "{} of data ({}) doesn't match one from file ({})."
            if dmin != field.min():
                dmin = field.min()
            if dmax != field.max():
                dmax = field.max()
#            assert np.allclose(dmin, field.min()), err_msg.format('Min', dmin,
#                                                                  field.min())
#            assert np.allclose(dmax, field.max()), err_msg.format('Max', dmax,
#                                                                  field.max())
            self.xmin = xmin
            self.ymin = ymin
            self.xdim = (xmax-xmin)/(nx-1)
            self.ydim = (ymax-ymin)/(ny-1)
            self.dmin = dmin
            self.dmax = dmax
            self.cols = ny
            self.rows = nx
            self.nullvalue = 1.701410009187828e+38
            self.data = np.ma.masked_equal(field.reshape(nx,ny), self.nullvalue)
        #x, y = gridder.regular(area, shape)
        #data = dict(file=fname, shape=shape, area=area, data=field, x=x, y=y)
        #return data

    def export_surfer(self, fname, flag = True ,file_format='binary'):
        """
        Export a surfer grid

        Parameters
        ----------
        fname : filename of grid dataset to export
        flag  : True  -  Output Grid Grid
                False -  Output Bak Grid Grid
        file_format : binary/b - output binary format
                      ascii/a - output ascii format
        """
        if file_format == 'binary' or file_format == 'b':
            self._export_surfer_binary(fname,flag)
        elif file_format == 'ascii' or file_format == 'a':
            self._export_surfer_ascii(fname,flag)

    def _export_surfer_ascii(self, fname, flag = True):
        """
        Export a surfer binary grid

        Parameters
        ----------
        fname : filename of grid dataset to export
        flag  : True  -  Output Grid Grid
                False -  Output Bak Grid Grid
        """
        xmax = self.xmin + (self.cols-1)*self.xdim
        ymax = self.ymin + (self.rows-1)*self.ydim
        with open(fname,'w') as fno:
            fno.write('DSAA\n')
            fno.write('{} {}\n'.format(self.cols,self.rows))
            fno.write('{} {}\n'.format(self.xmin,self.xmax))
            fno.write('{} {}\n'.format(self.ymin,self.ymax))
            if flag:
                fno.write('{} {}\n'.format(np.min(self.data),
                                           np.max(self.data))
                          )
                ntmp = 1.701410009187828e+38
                tmp = self.data.astype('f')
                tmp = tmp.filled(ntmp)
            else:
                fno.write('{} {}\n'.format(np.min(self.data0),
                                           np.max(self.data0))
                          )
                tmp = self.data0.astype('f')
            np.savetxt(fno,tmp)

    def _export_surfer_binary(self, fname, flag = True):
        """
        Export a surfer binary grid

        Parameters
        ----------
        fname : filename of grid dataset to export
        flag  : True  -  Output Grid Grid
                False -  Output Bak Grid Grid
        """
        fno = open(fname, 'wb')
        xmax = self.xmin + (self.cols-1)*self.xdim
        ymax = self.ymin + (self.rows-1)*self.ydim
        if flag:
            bintmp = struct.pack('cccchhdddddd', b'D', b'S', b'B', b'B',
                             self.cols, self.rows,
                             self.xmin, xmax,
                             self.ymin, ymax,
                             np.min(self.data),
                             np.max(self.data))
            fno.write(bintmp)
            ntmp = 1.701410009187828e+38
            tmp = self.data.astype('f')
            tmp = tmp.filled(ntmp)
        else:
            bintmp = struct.pack('cccchhdddddd', b'D', b'S', b'B', b'B',
                             self.cols, self.rows,
                             self.xmin, xmax,
                             self.ymin, ymax,
                             np.min(self.data0),
                             np.max(self.data0))
            fno.write(bintmp)
            tmp = self.data0.astype('f')
        #tmp = tmp[::-1]
        fno.write(tmp.tostring())
        fno.close()


    def export_ascii(self, fname):
        """
        Export Ascii file

        Parameters
        ----------
        data : grid Data
            dataset to export
        """
        fno = open(fname, 'w')

        fno.write("ncols \t\t\t" + str(self.cols))
        fno.write("\nnrows \t\t\t" + str(self.rows))
        fno.write("\nxllcorner \t\t\t" + str(self.xmin))
        fno.write("\nyllcorner \t\t\t" + str(self.ymin))
        fno.write("\ncellsize \t\t\t" + str(self.xdim))
        fno.write("\nnodata_value \t\t" + str(self.nullvalue))

        tmp = self.data.filled(self.nullvalue)

        for j in range(self.rows):
            fno.write("\n")
            for i in range(self.cols):
                fno.write(str(tmp[j, i]) + " ")

        fno.close()

if __name__ == "__main__":

 # 使用方法示例，masked numpy ndarray
 grd1=grddata()
 grd1.load_surfer(r'D:\demo\demogrid.grd')
 if np.ma.is_masked(grd1.data):
   grd1.fill_nulls()
   plt.imshow(grd1.data0)
 else:
   print('not null region in dataset')

 #d1=grd1.data
 #grd1.data=d1*d1
# v = d1.reshape(grd1.rows*grd1.cols)
# #gridder.interpolation.fill_nans(x, y, v, xp, yp, vp):
# plt.imshow(grd1.data)     #显示绘图结果
 grd1.export_surfer(r'D:\demo\demogrid3-blk.grd', flag = False)


