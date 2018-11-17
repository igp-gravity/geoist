#!/usr/bin/env python

# python 3 compatibility
from __future__ import print_function

# stdlib imports
import abc
import textwrap
import sys
import re

# third party imports
from .gridbase import Grid
from .dataset import DataSetException
from .geodict import GeoDict

import numpy as np
from scipy import interpolate
import shapely
from affine import Affine
from rasterio import features
from rasterio.warp import reproject, calculate_default_transform
from rasterio.crs import CRS
from rasterio.enums import Resampling
from osgeo import osr


def _center_projection(projection):
    parts = projection.split('+')[1:]
    newparts = []
    for part in parts:
        if 'lon_0' in part:
            part = 'lon_0=0 '
        newparts.append(part)
    newprojection = '+' + '+'.join(newparts)
    return newprojection


class Grid2D(Grid):
    reqfields = set(['xmin', 'xmax', 'ymin', 'ymax', 'dx', 'dy', 'nx', 'ny'])

    def __init__(self, data=None, geodict=None):
        """
        Construct a Grid object.

        :param data:
            A 2D numpy array (can be None).
        :param geodict:
            A GeoDict Object (or None) containing the following fields:
        :returns:
            A Grid2D object.
        :raises DataSetException:
           When data is not 2D or the number of rows and columns do not
           match the geodict.
        """
        if data is not None and geodict is not None:
            # complain if data is not 2D (squeeze 1d dimensions out)
            dims = data.shape
            if len(dims) != 2:
                raise DataSetException('Grid data must be 2D. Input data has '
                                       'shape of %s' % str(data.shape))
            ny, nx = dims
            if ny != geodict.ny or nx != geodict.nx:
                raise DataSetException('Input geodict does not match shape '
                                       'of input data.')
            self._geodict = geodict.copy()
            self._data = data.copy()
        else:
            self._data = None
            self._geodict = None

    @staticmethod
    def checkFirstColumnDuplicated(geodict):
        """Check to see if the first column in a file described by geodict
        is duplicated by the last column.

        :param geodict:
          GeoDict object which may have duplicate column.
        :returns:
          Tuple containing:
           - GeoDict object representing a grid with the last column removed.
           - Boolean indicating whether the last column was a duplicate.
        """
        first_column_duplicated = False
        cols_per_degree = 1/geodict.dx
        is_even = np.isclose(cols_per_degree, np.floor(cols_per_degree))
        is_global = (geodict.xmax - geodict.xmin) >= 360 or \
            (geodict.xmax-geodict.xmin) < 0
        newgeodict = geodict.copy()
        if is_even and is_global:
            if np.floor(cols_per_degree)*360 == geodict.nx-1:
                first_column_duplicated = True
                nx = geodict.nx-1
                xmax = geodict.xmin + (nx-1)*geodict.dx
                newgeodict = GeoDict({'xmin': geodict.xmin,
                                      'xmax': xmax,
                                      'ymin': geodict.ymin,
                                      'ymax': geodict.ymax,
                                      'nx': nx,
                                      'ny': geodict.ny,
                                      'dx': geodict.dx,
                                      'dy': geodict.dy})

        return (newgeodict, first_column_duplicated)

    @staticmethod
    def getDataRange(fgeodict, sampledict, first_column_duplicated=False,
                     padding=None):
        """For a given set of input bounds and information about a file,
        determine the rows and columns for bounds.

        :param fgeodict:
          GeoDict object for a given file.
        :param sampledict:
          Sampling GeoDict object.
        :param first_column_duplicated:
          Boolean indicating whether the last column in a file is a
          duplicate of the first column.
        :returns:
          Dictionary containing fields:
            - iulx1 Upper left X of first (perhaps only) segment.
            - iuly1 Upper left Y of first (perhaps only) segment.
            - ilrx1 Lower right X of first (perhaps only) segment.
            - ilry1 Lower right Y of first (perhaps only) segment.
            (if bounds cross 180 meridian...)
            - iulx2 Upper left X of second segment.
            - iuly2 Upper left Y of second segment.
            - ilrx2 Lower right X of second segment.
            - ilry2 Lower right Y of second segment.
        """
        data_range = {}

        # get the file values
        gxmin = fgeodict.xmin
        gxmax = fgeodict.xmax
        nx = fgeodict.nx
        ny = fgeodict.ny

        # make sure fgeodict and bounds are in the same longitude
        # system (-180 to 180, -360 to 0, or 0 to 360)
        xmin, xmax, ymax = (sampledict.xmin, sampledict.xmax,
                            sampledict.ymax)
        file_normal = gxmin >= -180 and gxmin <= 180 and \
            gxmax >= -180 and gxmax <= 180 and gxmax > gxmin

        bounds_negative = xmin < -180
        # bounds_360 = xmin >= 0 and xmin <= 360 and xmax >= 0 and
        # xmax <= 360 and xmax > gxmin
        bounds_360 = xmin >= 180 or xmax >= 180

        # if the file coordinates are between -180 and 180, inclusive...
        if file_normal:
            # if the input sampling bound xmin is less than -180...
            if bounds_negative:
                xmin += 360
                if xmax < -180:
                    xmax += 360
            if bounds_360:
                if xmin > 180:
                    xmin -= 360
                if xmax > 180:
                    xmax -= 360

        # these three values are the same regardless of whether we
        # cross meridian
        iuly1, iulx1 = fgeodict.getRowCol(ymax, xmin)
        ilrx1 = iulx1 + sampledict.nx
        ilry1 = iuly1 + sampledict.ny

        iulx2 = None
        ilrx2 = None
        iuly2 = None
        ilry2 = None
        if ilrx1 > nx:
            ilrx2 = ilrx1-nx
            ilrx1 = nx
            iulx2 = 0
            iuly2 = iuly1
            ilry2 = ilry1
        if ilry1 > ny:
            ilry1 = ny

        data_range['iulx1'] = int(iulx1)
        data_range['ilrx1'] = int(ilrx1)
        data_range['iuly1'] = int(iuly1)
        data_range['ilry1'] = int(ilry1)
        if iulx2 is not None:
            data_range['iulx2'] = int(iulx2)
            data_range['ilrx2'] = int(ilrx2)
            data_range['iuly2'] = int(iuly2)
            data_range['ilry2'] = int(ilry2)

        return data_range

    @staticmethod
    def verifyBounds(filegeodict, samplegeodict, resample=False):
        """Ensure that the two grids represented at least 1) intersect
        and 2) are aligned if resampling is True.

        :param filegeodict:
          GeoDict object describing file.
        :param samplegeodict:
          GeoDict object describing grid to use for sampling.
        :param resample:
          Boolean indicating whether we want to resample.
        :raises:
          DataSetException when geodicts do not intersect or if the grids
          are not aligned.
        """
        if samplegeodict is not None \
                and not filegeodict.intersects(samplegeodict):
            msg = 'Input samplegeodict must at least intersect with the '\
                  'file bounds'
            raise DataSetException(msg)
        # verify that if not resampling, the dimensions of the sampling
        # geodict must match the file.
        if resample is False and samplegeodict is not None:
            ddx = np.abs(filegeodict.dx - samplegeodict.dx)
            ddy = np.abs(filegeodict.dy - samplegeodict.dy)
            if ddx > GeoDict.EPS or ddy > GeoDict.EPS:
                raise DataSetException('File dimensions are different '
                                       'from sampledict dimensions.')

    @staticmethod
    def bufferBounds(samplegeodict, filegeodict, resample=False,
                     buffer_pixels=1, doPadding=False):
        """Buffer requested bounds out by buffer_pixels pixels, or edge
        of grid.

        Buffer pixels shoud be at filegeodict resolution.

        :param samplegeodict:
          GeoDict object describing grid to use for sampling.
        :param filegeodict:
          GeoDict object describing file.
        :param resample:
          Boolean indicating whether we want to resample.
        :param buffer_pixels:
          Number of pixels to buffer bounds in any possible direction.
        :returns:
          GeoDict which has been buffered by the appropriate number of pixels.
        """
        if not resample:
            return samplegeodict

        dx = filegeodict.dx
        dy = filegeodict.dy
        fxmin, fxmax = filegeodict.xmin, filegeodict.xmax
        fymin, fymax = filegeodict.ymin, filegeodict.ymax
        sxmin, sxmax = samplegeodict.xmin, samplegeodict.xmax
        symin, symax = samplegeodict.ymin, samplegeodict.ymax

        if not filegeodict.intersects(samplegeodict):
            if not doPadding:
                msg = 'Cannot buffer bounds when sampling grid is '\
                      'completely outside file grid, unless doPadding=True.'
                raise DataSetException(msg)
            else:
                return samplegeodict

        buffer_geo_x = buffer_pixels * dx
        buffer_geo_y = buffer_pixels * dy

        xmin = sxmin - buffer_geo_x
        xmax = sxmax + buffer_geo_x
        ymin = symin - buffer_geo_y
        ymax = symax + buffer_geo_y

        # assign the "most minimum" x value, taking 180 meridian into account
        is_180 = fxmin >= -180 and fxmax <= 180
        if xmin < fxmin:
            if is_180:
                if xmin > -180:
                    xmin = fxmin

        # assign the "most maximum" x value, taking 180 meridian into account
        if xmax > fxmax:
            if is_180:
                if xmax < 180:
                    xmax = fxmax

        if ymin < fymin:
            ymin = fymin

        if ymax > fymax:
            ymax = fymax

        # make sure that boundaries are on filegeodict boundaries -
        # if not, go outwards until we hit one
        ddxmin = xmin/dx
        if int(ddxmin) != ddxmin:
            xmin = np.floor(ddxmin) * dx
        ddxmax = xmax/dx
        if int(ddxmax) != ddxmax:
            xmax = np.ceil(ddxmax) * dx
        ddymin = ymin/dy
        if int(ddymin) != ddymin:
            ymin = np.floor(ddymin) * dy
        ddymax = ymax/dy
        if int(ddymax) != ddymax:
            ymax = np.ceil(ddymax) * dy

        geodict = GeoDict.createDictFromBox(xmin, xmax, ymin, ymax, dx, dy,
                                            inside=True)
        return geodict

    @staticmethod
    def padGrid(data, geodict, pad_dict):
        """Pad input data array with pixels specified by pad_dict on each side.

        :param data:
          2D numpy array of data.
        :param geodict:
          GeoDict object describing data.
        :param pad_dict:
          A dictionary containing fields:
            - padleft The number of padding pixels on the left edge.
            - padright The number of padding pixels on the right edge.
            - padbottom The number of padding pixels on the bottom edge.
            - padtop The number of padding pixels on the top edge.
        :returns:
          Tuple of (data,geodict) where data has been padded and geodict
          represents new padded data.
        """
        if pad_dict['padleft'] == 0 and pad_dict['padright'] == 0 and \
                pad_dict['padbottom'] == 0 and pad_dict['padtop'] == 0:
            return (data, geodict)
        newdata = data.copy()
        ny, nx = newdata.shape
        dx, dy = geodict.dx, geodict.dy
        # we're padding with inf so that we don't interfere with other nan
        # pixels in the data
        leftcols = np.ones((ny, pad_dict['padleft'])) * np.inf
        newdata = np.hstack((leftcols, newdata))
        ny, nx = newdata.shape
        rightcols = np.ones((ny, pad_dict['padright'])) * np.inf
        newdata = np.hstack((newdata, rightcols))
        ny, nx = newdata.shape
        bottomrows = np.ones((pad_dict['padbottom'], nx)) * np.inf
        newdata = np.vstack((newdata, bottomrows))
        ny, nx = newdata.shape
        toprows = np.ones((pad_dict['padtop'], nx)) * np.inf
        newdata = np.vstack((toprows, newdata))

        ny, nx = newdata.shape

        # get the shapes of the pads - sometimes these can have a shape
        # of (3,0) and a length of 3
        leftheight, leftwidth = leftcols.shape
        rightheight, rightwidth = rightcols.shape
        bottomheight, bottomwidth = bottomrows.shape
        topheight, topwidth = toprows.shape

        newxmin = geodict.xmin - leftwidth*dx
        newxmax = geodict.xmax + rightwidth*dx
        newymin = geodict.ymin - bottomheight*dy
        newymax = geodict.ymax + topheight*dy
        newdict = GeoDict({'xmin': newxmin,
                           'xmax': newxmax,
                           'ymin': newymin,
                           'ymax': newymax,
                           'nx': nx,
                           'ny': ny,
                           'dx': dx,
                           'dy': dy}, adjust='res')
        return (newdata, newdict)

    @staticmethod
    def getPadding(filegeodict, samplegeodict, doPadding=False):
        """Determine how many pixels of padding there need to be on each
        side of requested grid.

        :param filegeodict:
          GeoDict object specifying the spatial information from a source file.
        :param samplegeodict:
          GeoDict object specifying the spatial information for a desired
          sampling regime.
        :param resampling:
          Boolean indicating that user wants to resample the data from
          the file to the samplegeodict.
        :raises DataSetException:
          When resampling is False and filegeodict and samplegeodict are
          not pixel aligned.
        :returns:
          A dictionary containing fields:
            - padleft The number of padding pixels on the left edge.
            - padright The number of padding pixels on the right edge.
            - padbottom The number of padding pixels on the bottom edge.
            - padtop The number of padding pixels on the top edge.
        """
        if not doPadding:
            pad_dict = {'padleft': 0,
                        'padright': 0,
                        'padbottom': 0,
                        'padtop': 0}
        else:
            # get pad left columns - go outside specified bounds if not
            # exact edge
            pxmin, pxmax, pymin, pymax = (samplegeodict.xmin,
                                          samplegeodict.xmax,
                                          samplegeodict.ymin,
                                          samplegeodict.ymax)
            gxmin, gxmax, gymin, gymax = (filegeodict.xmin,
                                          filegeodict.xmax,
                                          filegeodict.ymin,
                                          filegeodict.ymax)
            dx, dy = (filegeodict.dx, filegeodict.dy)

            padleftcols = int(np.ceil((gxmin - pxmin)/dx))
            padrightcols = int(np.ceil((pxmax - gxmax)/dx))
            padbottomrows = int(np.ceil((gymin - pymin)/dy))
            padtoprows = int(np.ceil((pymax - gymax)/dy))

            # if any of these are negative, set them to zero
            if padleftcols < 0:
                padleftcols = 0
            if padrightcols < 0:
                padrightcols = 0
            if padbottomrows < 0:
                padbottomrows = 0
            if padtoprows < 0:
                padtoprows = 0
            pad_dict = {'padleft': padleftcols,
                        'padright': padrightcols,
                        'padbottom': padbottomrows,
                        'padtop': padtoprows}
        return pad_dict

    def __repr__(self):
        """
        String representation of a Grid2D object.
        :returns:
          String containing description of Grid2D object.
        """
        xmin, xmax, ymin, ymax = (self._geodict.xmin, self._geodict.xmax,
                                  self._geodict.ymin, self._geodict.ymax)
        ny, nx = self._data.shape
        dx, dy = (self._geodict.dx, self._geodict.dy)
        zmin = np.nanmin(self._data)
        zmax = np.nanmax(self._data)
        rstr = '''<%s Object:
        ny: %i
        nx: %i
        xmin: %.4f
        xmax: %.4f
        ymin: %.4f
        ymax: %.4f
        dx: %.4f
        dy: %.4f
        zmin: %.6f
        zmax: %.6f>''' % (self.__class__.__name__, ny, nx, xmin, xmax,
                          ymin, ymax, dx, dy, zmin, zmax)
        parts = rstr.split('\n')
        newrstr = '\n'.join([p.strip() for p in parts])
        return textwrap.dedent(newrstr)

    @classmethod
    def _createSampleData(self, M, N):
        """Used for internal testing - create an NxN grid with lower left
        corner at 0.5,0.5, dx/dy = 1.0.

        :param M:
           Number of rows in output grid
        :param N:
           Number of columns in output grid
        :returns:
           GMTGrid object where data values are an NxN array of values
           from 0 to N-squared minus 1, and geodict
           lower left corner is at 0.5/0.5 and cell dimensions are 1.0.
        """
        data = np.arange(0, M*N).reshape(M, N)
        # arange gives int64 by default, not supported by netcdf3
        data = data.astype(np.int32)
        xvar = np.arange(0.5, 0.5+N, 1.0)
        yvar = np.arange(0.5, 0.5+M, 1.0)
        geodict = {'ny': M,
                   'nx': N,
                   'xmin': 0.5,
                   'xmax': xvar[-1],
                   'ymin': 0.5,
                   'ymax': yvar[-1],
                   'dx': 1.0,
                   'dy': 1.0}
        gd = GeoDict(geodict)
        return (data, gd)

    @abc.abstractmethod
    def getFileGeoDict(filename):
        # this should return a geodict and a boolean indicating whether
        # the first
        # column is duplicated at the end of each row (some data providers
        # do this).
        raise NotImplementedError('Load method not implemented in base class')

    @abc.abstractmethod
    def readFile(filename, data_range):
        """Read in data from the given file, at the pixels specified in
        data_range.

        :param filename:
          Name of file to read.
        :param data_range:
          Dictionary containing fields:
            - iulx1 Upper left X of first (perhaps only) segment.
            - iuly1 Upper left Y of first (perhaps only) segment.
            - ilrx1 Lower right X of first (perhaps only) segment.
            - ilry1 Lower right Y of first (perhaps only) segment.
            (if bounds cross 180 meridian...)
            - iulx2 Upper left X of second segment.
            - iuly2 Upper left Y of second segment.
            - ilrx2 Lower right X of second segment.
            - ilry2 Lower right Y of second segment.
        """
        raise NotImplementedError('readFile is only implemented in Grid2D '
                                  'subclasses.')

    @classmethod
    def copyFromGrid(cls, grid):
        """
        Copy constructor - can be used to create an instance of any
        Grid2D subclass from another.

        :param grid:
          Any Grid2D instance.
        :returns:
          A copy of the data in that input Grid2D instance.
        """
        if not isinstance(grid, Grid2D):
            raise DataSetException('Input to copyFromGrid must be an '
                                   'instance of a Grid2D object (inc. '
                                   'subclasses)')
        return cls(grid.getData(), grid.getGeoDict())

    # This should be a @classmethod in subclasses
    @abc.abstractmethod
    def save(self, filename):
        # would we ever want to save a subset of the data?
        raise NotImplementedError('Save method not implemented in base class')

    @classmethod
    def _createSections(self, bounds, geodict, firstColumnDuplicated,
                        isScanLine=False):
        """Given a grid that goes from -180 to 180 degrees, figure out
        the two pixel regions that up both sides of the subset.

        :param bounds:
           Tuple of (xmin,xmax,ymin,ymax)
        :param geodict:
           Geodict dictionary
        :param firstColumnDuplicated:
          Boolean indicating whether this is a global data set where the
          first and last columns are identical.
        :param isScanLine:
          Boolean indicating whether this array is in scan line order
          (pixel[0,0] is the geographic upper left).
        :returns:
          Two tuples of 4 elements each - (iulx,iuly,ilrx,ilry). The
          first tuple defines the pixel offsets for the left
          side of the subsetted region, and the second tuple
          defines the pixel offsets for the right side.
        """
        (bxmin, bxmax, bymin, bymax) = bounds
        ulx = geodict.xmin
        uly = geodict.ymax
        dx = geodict.dx
        dy = geodict.dy
        nx = geodict.nx
        # section 1
        iulx1 = int(np.floor((bxmin - ulx)/dx))
        ilrx1 = int(nx)

        if not isScanLine:
            iuly1 = int(np.ceil((uly - bymax)/dy))
            ilry1 = int(np.floor((uly - bymin)/dy)) + 1
        else:
            ilry1 = int(np.ceil((uly - bymin)/dy))
            iuly1 = int(np.floor((uly - bymax)/dy)) + 1

        # section 2
        iulx2 = 0
        ilrx2 = int(np.ceil((bxmax - ulx)/dx)) + 1
        iuly2 = iuly1
        ilry2 = ilry1

        if firstColumnDuplicated:
            ilrx1 -= 1

        region1 = (iulx1, iuly1, ilrx1, ilry1)
        region2 = (iulx2, iuly2, ilrx2, ilry2)
        return(region1, region2)

    def getData(self):
        """Return a reference to the data inside the Grid.

        :returns:
          A reference to a 2D numpy array.
        """
        return self._data  # should we return a copy of the data?

    def setData(self, data):
        """Set the data inside the Grid.

        :param data:
          A 2D numpy array.
        :raises:
          DataSetException if the number of rows and columns do not match
          the the internal GeoDict, or if the input
          is not a numpy array.
        """
        if not isinstance(data, np.ndarray):
            raise DataSetException('setData() input is not a numpy array.')

        if len(data.shape) != 2:
            raise DataSetException('setData() input array must have '
                                   'two dimensions.')

        m, n = data.shape
        if m != self._geodict.ny or n != self._geodict.nx:
            raise DataSetException('setData() input array must match rows '
                                   'and columns of existing data.')
        self._data = data

    def getGeoDict(self):
        """
        Return a reference to the geodict inside the Grid.
        :returns:
          A reference to a dictionary (see constructor).
        """
        return self._geodict.copy()

    def getBounds(self):
        """
        Return the lon/lat range of the data.

        :returns:
           Tuple of (lonmin,lonmax,latmin,latmax)
        """
        return (self._geodict.xmin, self._geodict.xmax, self._geodict.ymin,
                self._geodict.ymax)

    def applyNaN(self, force=False):
        """Apply no data value to internal data, cast to float if necessary.

        Intelligently cast data in grid to be able to handle NaN values.

        Usage:
        Integer data with a precision of 16 bits or less
        will be cast to 32 bit floating point.

        Integer data with precision of 32 bits will be cast to 32 bit
        floating point if maximum/minimum values can be cast without
        losing precision.

        Integer data with precision of 64 bits or greater will be cast to
        64 bit floating point if maximum/minimum values can be cast without
        losing precision. Otherwise, this method will raise an OverflowError
        unless the force option is set to True.

        :param force:
          Boolean indicating whether to override OverflowError (see Usage).
        """
        nodata = self._geodict.nodata
        if nodata is None or np.isnan(nodata) or np.isnan(self._data).any():
            return
        isint = 'int' in str(self._data.dtype)
        precision = self._data.dtype.itemsize
        if not isint:
            self._data[self._data == nodata] = np.nan
        if isint:
            if precision <= 2:
                self._data = self._data.astype(np.float32)
            elif precision <= 4:
                dmax = self._data.max()
                dmin = self._data.min()
                fdmax = np.float32(dmax)
                fdmin = np.float32(dmin)
                if dmax == fdmax and dmin == fdmin:
                    self._data = self._data.astype(np.float32)
                else:
                    self._data = self._data.astype(np.float64)
            else:
                dmax = self._data.max()
                dmin = self._data.min()
                fdmax = np.float64(dmax)
                fdmin = np.float64(dmin)
                if dmax == fdmax and dmin == fdmin:
                    if not force:
                        raise OverflowError(
                            'Data cannot be represented as 64-bit float.')
                self._data = self._data.astype(np.float64)
            # if we've gotten this far, we have upcasted successfully
            self._data[self._data == nodata] = np.nan

    def subdivide(self, finerdict, cellFill='max'):
        """Subdivide the cells of the host grid into finer-resolution cells.

        :param finerdict:
          GeoDict object defining a grid with a finer resolution than the
          host grid.
        :param cellFill:
          String defining how to fill cells that span more than one host
          grid cell.
          Choices are:
            'max': Choose maximum value of host grid cells.
            'min': Choose minimum value of host grid cells.
            'mean': Choose mean value of host grid cells.
        :returns:
          Grid2D instance with host grid values subdivided onto finer grid.
        :raises DataSetException:
          When finerdict is not a) finer resolution or b) does not
          intersect.x or cellFill is not valid.
        """
        fillvals = ['min', 'max', 'mean']
        if cellFill not in fillvals:
            raise DataSetException('cellFill input must be one of %s.' %
                                   fillvals)
        if finerdict.dx >= self._geodict.dx or finerdict.dy >= self._geodict.dy:
            raise DataSetException('subdivide() input GeoDict must be finer '
                                   'resolution than host grid.')
        if not finerdict.intersects(self._geodict):
            raise DataSetException('subdivide() input GeoDict must intersect '
                                   'host grid.')

        # things are simple if the host grid cell dx/dy are a multiple
        # of finer grid dx/dy and are
        # aligned in the sense that every host grid cell edge matches
        # an edge of finer grid cell.
        resXMultiple = self._geodict.dx/finerdict.dx == \
            int(self._geodict.dx/finerdict.dx)
        resYMultiple = self._geodict.dy/finerdict.dy == \
            int(self._geodict.dy/finerdict.dy)
        # this stuff below may not be right...?
        dxmin = (self._geodict.xmin-finerdict.xmin)/finerdict.dx
        isXAligned = np.isclose(dxmin, int(dxmin))
        dymin = (self._geodict.ymin-finerdict.ymin)/finerdict.dy
        isYAligned = np.isclose(dymin, int(dymin))
        isAligned = resXMultiple and resYMultiple and isXAligned and isYAligned
        finedata = np.ones((finerdict.ny, finerdict.nx),
                           dtype=self._data.dtype)*np.nan
        if isAligned:
            for i in range(0, self._geodict.ny):
                for j in range(0, self._geodict.nx):
                    cellvalue = self._data[i, j]
                    # what is the longitude of the first finer cell inside
                    # the host cell?
                    # coordinates of center of host cell
                    clat, clon = self.getLatLon(i, j)
                    # get the left edge of the cell
                    fleftlon = clon - (self._geodict.dx/2) + finerdict.dx/2
                    ftoplat = clat + (self._geodict.dy/2) - finerdict.dy/2
                    frightlon = clon + (self._geodict.dx/2) - finerdict.dx/2
                    fbottomlat = clat - (self._geodict.dy/2) + finerdict.dy/2
                    itop, jleft = finerdict.getRowCol(ftoplat, fleftlon)
                    itop = itop
                    jleft = jleft
                    ibottom, jright = finerdict.getRowCol(fbottomlat,
                                                          frightlon)
                    ibottom = ibottom
                    jright = jright
                    finedata[itop:ibottom+1, jleft:jright+1] = cellvalue
        else:
            for i in range(0, self._geodict.ny):
                for j in range(0, self._geodict.nx):
                    cellvalue = self._data[i, j]
                    # get the indices of all cells that are
                    # completely contained inside this one.
                    # coordinates of center of host cell
                    clat, clon = self.getLatLon(i, j)
                    # what is the longitude of of our first approximated
                    # left edge fine
                    # cell that is contained by host cell?
                    fleftlon = clon - self._geodict.dx/2.0 + finerdict.dx/2
                    frightlon = clon + self._geodict.dx/2.0 - finerdict.dx/2
                    jleft = int(np.ceil((fleftlon - finerdict.xmin) /
                                        finerdict.dx))
                    jright = int(np.floor((frightlon - finerdict.xmin) /
                                          finerdict.dx))

                    # what is the latitude of of our first approximated
                    # bottom edge fine
                    # cell that is contained by host cell?
                    fbottomlat = clat - self._geodict.dy/2.0 + finerdict.dy/2
                    ftoplat = clat + self._geodict.dy/2.0 - finerdict.dy/2
                    ibottom = int(np.floor((finerdict.ymax - fbottomlat) /
                                           finerdict.dy))
                    itop = int(np.ceil((finerdict.ymax - ftoplat) /
                                       finerdict.dy))
                    # ibottom = int(np.ceil((fbottomlat - finerdict.ymin) /
                    # finerdict.dy))
                    # itop = int(np.floor((ftoplat - finerdict.ymin) /
                    # finerdict.dy))

                    finedata[itop:ibottom+1, jleft:jright+1] = cellvalue
                    # now what do I do about cells that aren't completely
                    # contained?

            # we have to now find all rows/columns where there are NaN
            # values and deal with them
            # accordingly - let's look at output rows first, looking for a
            # row that is all NaN
            # and doesn't have an all NaN row above or below it.

            # np.minimum and maximum throw warnings with NaN values, even
            # the behavior with respect
            # to those is clearly documented.  So, let's ignore those
            # warnings in this
            # section where we use those methods.
            with np.errstate(invalid='ignore'):
                colidx = finerdict.nx//2
                while colidx > -1:
                    col = finedata[:, colidx]
                    if not np.isnan(col).all():
                        nanrows = np.where(np.isnan(col))
                        break
                    colidx -= 1
                for i in nanrows[0]:
                    if i == 0 or i == finerdict.ny-1:
                        continue
                    if cellFill == 'min':
                        finedata[i, :] = np.minimum(finedata[i-1, :],
                                                    finedata[i+1, :])
                    elif cellFill == 'max':
                        finedata[i, :] = np.maximum(finedata[i-1, :],
                                                    finedata[i+1, :])
                    else:  # cellFill == 'mean':
                        finedata[i, :] = (finedata[i-1, :] +
                                          finedata[i+1, :]) / 2.0
                # now look at output columns
                rowidx = finerdict.ny//2
                while rowidx > -1:
                    row = finedata[rowidx, :]
                    if not np.isnan(row).all():
                        nancols = np.where(np.isnan(row))
                        break
                    rowidx -= 1
                for j in nancols[0]:
                    if j == 0 or j == finerdict.nx-1:
                        continue
                    if cellFill == 'min':
                        finedata[:, j] = np.minimum(finedata[:, j-1],
                                                    finedata[:, j+1])
                    elif cellFill == 'max':
                        finedata[:, j] = np.maximum(finedata[:, j-1],
                                                    finedata[:, j+1])
                    else:  # cellFill == 'mean':
                        finedata[:, j] = (finedata[:, j-1] +
                                          finedata[:, j+1]) / 2.0

        finegrid = Grid2D(finedata, finerdict)
        return finegrid

    def cut(self, xmin, xmax, ymin, ymax, align=False):
        """Cut out a section of Grid and return it.

        :param xmin: Longitude coordinate of upper left pixel, must be
            aligned with Grid.
        :param xmax: Longitude coordinate of lower right pixel, must be
            aligned with Grid.
        :param ymin: Latitude coordinate of upper left pixel, must be
            aligned with Grid.
        :param ymax: Latitude coordinate of lower right pixel, must be
            aligned with Grid.
        :param align: Boolean indicating whether input boundaries
            should be modified to align with host grid.
        """
        td1 = GeoDict.createDictFromBox(xmin, xmax, ymin, ymax,
                                        self._geodict.dx,
                                        self._geodict.dy, inside=True)
        td = None
        if not td1.isAligned(self._geodict):
            if not align:
                raise DataSetException('Input bounds must align with this '
                                       'grid.')
            else:
                td = self._geodict.getAligned(td1, inside=True)
        else:
            td = td1.copy()
        if not self._geodict.contains(td):
            raise DataSetException('Input bounds must be completely '
                                   'contained by this grid.')
        uly, ulx = self._geodict.getRowCol(td.ymax, td.xmin)
        lry, lrx = self._geodict.getRowCol(td.ymin, td.xmax)
        data = self._data[uly:lry+1, ulx:lrx+1]
        grid = Grid2D(data, td)
        return grid

    def getValue(self, lat, lon, method='nearest', default=None):
        # return nearest neighbor value
        """Return numpy array at given latitude and longitude (using
        nearest neighbor).

        :param lat:
           Latitude (in decimal degrees) of desired data value.
        :param lon:
           Longitude (in decimal degrees) of desired data value.
        :param method:
           Interpolation method, one of ('nearest','linear','cubic','quintic')
        :param default:
           Default value to return when lat/lon is outside of grid bounds.
        :return:
           Value at input latitude,longitude position.
        :raises DataSetException:
          When lat/lon is outside of bounds and default is None.
        """

        if method == 'nearest':
            row, col = self.getRowCol(lat, lon)
        else:
            raise NotImplementedError('nearest is the only interpolation '
                                      'method currently supported.')
        ny, nx = self._data.shape
        if isinstance(row, np.ndarray):
            outidx = np.where((row < 0) | (row > ny-1) | (col < 0) |
                              (col > nx-1))
            inidx = np.where((row >= 0) & (row <= ny-1) & (col >= 0) &
                             (col <= nx-1))
            value = np.empty_like(row).astype(self._data.dtype)
            if len(outidx[0]):
                if default is None:
                    msg = 'One of more of your lat/lon values is outside '\
                        'Grid boundaries: %s' % (str(self.getBounds()))
                    raise DataSetException(msg)
                value[outidx] = default
            value[inidx] = self._data[row[inidx], col[inidx]]
        else:
            if (row < 0 or row > ny-1 or col < 0 or col > nx-1):
                if default is None:
                    msg = 'Your lat/lon value is outside Grid boundaries: %s'\
                        % (str(self.getBounds()))
                    raise DataSetException(msg)
                else:
                    return default
            value = self._data[row, col]
        return value

    def getLatLon(self, row, col):
        """Return geographic coordinates (lat/lon decimal degrees) for
        given data row and column.

        :param row:
           Row dimension index into internal data array.
        :param col:
           Column dimension index into internal data array.
        :returns:
           Tuple of latitude and longitude.
        """
        return self._geodict.getLatLon(row, col)

    def getRowCol(self, lat, lon, returnFloat=False):
        """Return data row and column from given geographic coordinates
        (lat/lon decimal degrees).

        :param lat:
           Input latitude.
        :param lon:
           Input longitude.
        :param returnFloat:
           Boolean indicating whether floating point row/col coordinates
           should be returned.
        :returns:
           Tuple of row and column.
        """
        return self._geodict.getRowCol(lat, lon, returnFloat)

    def _getInterpCoords(self, geodict):
        # translate geographic coordinates to 2 1-D arrays of X and Y
        # pixel coordinates
        # remember that pixel coordinates are (0,0) at the top left and
        # increase going down and to the right
        # geographic coordinates are (xmin,ymin) at the bottom left and
        # increase going up and to the right

        # we need to check if these two grids (host and sample) are in
        # the same longitude system
        # if not, adjust one to the other
        hostxmin = self._geodict.xmin
        hostxmax = self._geodict.xmax
        hostymin = self._geodict.ymin
        hostymax = self._geodict.ymax
        host_normal = hostxmin >= -180 and hostxmin <= 180 and \
            hostxmax >= -180 and hostxmax <= 180
        host_360 = hostxmin >= 180 or hostxmax >= 180

        samplexmin = geodict.xmin
        samplexmax = geodict.xmax
        sampleymin = geodict.ymin
        sampleymax = geodict.ymax
        sample_normal = samplexmin >= -180 and samplexmin <= 180 and \
            samplexmax >= -180 and samplexmax <= 180
        sample_360 = samplexmin >= 180 or samplexmax >= 180

        # get everything into negative mode
        if host_normal:
            if hostxmin > hostxmax:
                hostxmin -= 360
        if host_360:
            if hostxmin > 180:
                hostxmin -= 360
            if hostxmax > 180:
                hostxmax -= 360
        if sample_normal:
            if samplexmin > samplexmax:
                samplexmin -= 360
        if sample_360:
            if samplexmin > 180:
                samplexmin -= 360
            if samplexmax > 180:
                samplexmax -= 360

        # extract the geographic information about the grid we're sampling to
        sampleny = geodict.ny
        samplenx = geodict.nx

        hostdx = self._geodict.dx
        hostdy = self._geodict.dy

        # make sure that the grid we're resampling TO is completely
        # contained by host
        if samplexmin < hostxmin or samplexmax > hostxmax or \
                sampleymin < hostymin or sampleymax > hostymax:
            raise DataSetException('Grid you are resampling TO is not '
                                   'completely contained by base grid.')

        gxi = np.linspace(samplexmin, samplexmax, num=samplenx)
        gyi = np.linspace(sampleymin, sampleymax, num=sampleny)

        xi = (gxi - hostxmin)/hostdx
        yi = np.array(sorted(((hostymax - gyi)/hostdy)))

        return (xi, yi)

    def interpolate2(self, geodict, method='linear'):
        """
        Given a geodict specifying another grid extent and resolution,
        resample current grid to match.

        This method uses the rasterio reproject method instead of scipy.
        In tests with small grids
        that can be replicated easily by hand, the results from this
        method match interpolateToGrid.
        Limited tests with larger, random grids indicate some differences,
        the extent of which
        has not been documented.  These same limited tests indicate
        this method is 5 to 7 times
        faster than interpolateToGrid.

        :param geodict:
            geodict dictionary from another grid whose extents are
            inside the extent of this grid.
        :param method:
            Optional interpolation method - ['linear', 'cubic','nearest']
        :raises DataSetException:
           If the Grid object upon which this function is being called is
           not completely contained by the grid to which this Grid is being
           resampled.
        :raises DataSetException:
           If the method is not one of ['nearest','linear','cubic']
           If the resulting interpolated grid shape does not match input
           geodict.
        :returns:
          A new instance of the Grid2D class or subclass with interpolated
          data.
        """
        resampling = None
        if method == 'linear':
            resampling = Resampling.bilinear
        elif method == 'cubic':
            resampling = Resampling.cubic
        elif method == 'nearest':
            resampling = Resampling.nearest
        else:
            raise DataSetException('Unknown interpolation method %s' % method)

        destination = np.zeros((geodict.ny, geodict.nx))
        src_transform = Affine.from_gdal(self._geodict.xmin -
                                         self._geodict.dx/2.0,
                                         self._geodict.dx,
                                         0.0,  # x rotation, not used by us
                                         self._geodict.ymax +
                                         self._geodict.dy/2.0,
                                         0.0,  # y rotation, not used by us
                                         -1 * self._geodict.dy)  # their dy
        # is negative
        dst_transform = Affine.from_gdal(geodict.xmin - geodict.dx/2.0,
                                         geodict.dx,
                                         0.0,  # x rotation, not used by us
                                         geodict.ymax + geodict.dy/2.0,
                                         0.0,  # y rotation, not used by us
                                         -1*geodict.dy)  # their dy is negative
        src_crs = CRS().from_string(GeoDict.DEFAULT_PROJ4).to_dict()
        dst_crs = CRS().from_string(GeoDict.DEFAULT_PROJ4).to_dict()
        reproject(self._data.astype(np.float64), destination,
                  src_transform=src_transform,
                  src_crs=src_crs,
                  src_nodata=np.nan,
                  dst_transform=dst_transform,
                  dst_crs=dst_crs,
                  dst_nodata=np.nan,
                  resampling=resampling)

        return self.__class__(destination, geodict)

    def interpolateToGrid(self, geodict, method='linear'):
        """
        Given a geodict specifying another grid extent and resolution,
        resample current grid to match.

        :param geodict:
            geodict dictionary from another grid whose extents are inside
            the extent of this grid.
        :param method:
            Optional interpolation method - ['linear', 'cubic','nearest']
        :raises DataSetException:
           If the Grid object upon which this function is being called is
           not completely contained by the grid to which this Grid is being
           resampled.
        :raises DataSetException:
           If the method is not one of ['nearest','linear','cubic']
           If the resulting interpolated grid shape does not match input
           geodict.
        :returns:
          A new instance of the Grid2D class or subclass with interpolated
          data.
        """
        if method not in ['linear', 'cubic', 'nearest']:
            raise DataSetException('Resampling method must be one of '
                                   '"linear", "cubic","nearest"')
        ny, nx = (geodict.ny, geodict.nx)
        xi, yi = self._getInterpCoords(geodict)

        # now using scipy interpolate functions
        baserows, basecols = self._geodict.ny, self._geodict.nx
        basex = np.arange(0, basecols)  # base grid PIXEL coordinates
        basey = np.arange(0, baserows)
        newdata = None
        if method in ['linear', 'cubic']:
            hasnan = np.isnan(self._data).any()
            hasinf = np.isinf(self._data).any()
            if not hasnan and not hasinf:
                # at the time of this writing, interp2d does not support
                # NaN values at all.
                f = interpolate.interp2d(basex, basey, self._data, kind=method)
                # self._data = f(xi,yi)
                newdata = f(xi, yi)
                newdata = newdata.reshape((yi.size, xi.size))
            else:
                # is Nan values are present, use griddata (slower by ~2
                # orders of magnitude but supports NaN).
                xi, yi = np.meshgrid(xi, yi)
                newrows, newcols = xi.shape
                xi = xi.flatten()
                yi = yi.flatten()
                xnew = np.zeros((len(xi), 2))
                xnew[:, 0] = xi
                xnew[:, 1] = yi
                basex, basey = np.meshgrid(basex, basey)
                basex = basex.flatten()
                basey = basey.flatten()
                xold = np.zeros((len(basex), 2))
                xold[:, 0] = basex
                xold[:, 1] = basey
                # self._data = interpolate.griddata(xold,
                # self._data.flatten(),xnew,method=method)
                # self._data = self._data.reshape((newrows,newcols))
                newdata = interpolate.griddata(xold, self._data.flatten(),
                                               xnew, method=method)
                newdata = newdata.reshape((newrows, newcols))
        else:
            x, y = np.meshgrid(basex, basey)
            # in Python2, list doesn't do anything
            # in python3, it makes result of zip from iterator into list
            xy = list(zip(x.flatten(), y.flatten()))
            f = interpolate.NearestNDInterpolator(xy, self._data.flatten())
            newrows = geodict.ny
            newcols = geodict.nx
            xi = np.tile(xi, (newrows, 1))
            yi = np.tile(yi.reshape(newrows, 1), (1, newcols))
            # self._data = f(list(zip(xi.flatten(),yi.flatten())))
            # self._data = self._data.reshape(xi.shape)
            newdata = f(list(zip(xi.flatten(), yi.flatten())))
            newdata = newdata.reshape(xi.shape)

        ny, nx = geodict.ny, geodict.nx
        # dims = self._data.shape
        dims = newdata.shape
        ny_new = dims[0]
        nx_new = dims[1]
        if ny_new != ny or nx_new != nx:
            msg = "Interpolation failed!  Results (%i,%i) don't match "\
                "(%i,%i)!" % (ny_new, nx_new, ny, nx)
            raise DataSetException(msg)
        # now the extents and resolution of the two grids should be
        # identical...
        gdict = {'ny': geodict.ny,
                 'nx': geodict.nx,
                 'xmin': geodict.xmin,
                 'xmax': geodict.xmax,
                 'ymin': geodict.ymin,
                 'ymax': geodict.ymax,
                 'dx': geodict.dx,
                 'dy': geodict.dy}
        # self._geodict = GeoDict(gdict)
        newdict = GeoDict(gdict)
        return self.__class__(newdata, newdict)

    @classmethod
    def rasterizeFromGeometry(cls, shapes, geodict, burnValue=1.0,
                              fillValue=np.nan,
                              mustContainCenter=False, attribute=None):
        """
        Create a Grid2D object from vector shapes, where the presence
        of a shape
        (point, line, polygon) inside a cell turns that cell "on".

        :param shapes:
          One of:
            - One shapely geometry object (Point, Polygon, etc.) or a
              sequence of such objects
            - One GeoJSON like object or sequence of such objects.
              (http://geojson.org/)
            - A tuple of (geometry,value) or sequence of (geometry,value).
        :param geodict:
          GeoDict object which defines the grid onto which the shape values
          should be "burned".
        :param burnValue:
          Optional value which will be used to set the value of the pixels
          if there is no
          value in the geometry field.
        :param fillValue:
          Optional value which will be used to fill the cells not touched
          by any geometry.
        :param mustContainCenter:
          Optional boolean which indicates whether the geometry must touch
          the center of the cell or merely be inside the cell in order to
          set the value.
        :raises DataSetException:
          When geometry input is not a subclass of
          shapely.geometry.base.BaseGeometry.
        :returns:
          Grid2D object.
        This method is a thin wrapper around rasterio->features->rasterize(),
        documented here:
        https://github.com/mapbox/rasterio/blob/master/docs/features.rst

        which is itself a Python wrapper around the functionality found
        in gdal_rasterize, documented here:
        http://www.gdal.org/gdal_rasterize.html
        """
        # check the type of shapes
        # features.rasterize() documentation says this:
        # iterable of (geometry, value) pairs or iterable over
        # geometries. `geometry` can either be an object that implements
        # the geo interface or GeoJSON-like object.

        # create list of allowable types
        if sys.version_info.major == 2:
            types = (int, float, long)
        else:
            types = (int, float)

        # figure out whether this is a single shape or a sequence of shapes
        isOk = False
        isShape = False
        if isinstance(shapes, shapely.geometry.base.BaseGeometry):
            isOk = True
            isShape = True
        elif len(shapes) and isinstance(shapes[0],
                                        shapely.geometry.base.BaseGeometry):
            isOk = True
            isShape = True
        elif isinstance(shapes, dict) and 'geometry' in shapes and \
                'properties' in shapes:
            isOk = True
        elif len(shapes) and isinstance(shapes[0], dict) and \
                'geometry' in shapes[0] and 'properties' in shapes[0]:
            isOk = True
        else:
            pass
        if not isOk:
            raise DataSetException('shapes must be a single shapely object '
                                   'or sequence of them, or single Geo-JSON '
                                   'like-object')

        if not isShape:
            shapes2 = []
            for shape in shapes:
                geometry = shape['geometry']
                props = shape['properties']
                if attribute is not None:
                    if attribute not in props:
                        raise DataSetException('Input shapes do not have '
                                               'attribute "%s".' % attribute)
                    value = props[attribute]
                    if not isinstance(value, types):
                        raise DataSetException('value from input shapes '
                                               'object is not a number')
                else:
                    value = burnValue
                shapes2.append((geometry, value))
            shapes = shapes2

        xmin, xmax, ymin, ymax = (geodict.xmin, geodict.xmax,
                                  geodict.ymin, geodict.ymax)
        dx, dy = (geodict.dx, geodict.dy)

        if xmax < xmin:
            xmax += 360
        xvar = np.arange(xmin, xmax+(dx*0.1), dx)
        yvar = np.arange(ymin, ymax+(dy*0.1), dy)
        nx = len(xvar)
        ny = len(yvar)

        # the rasterize function assumes a pixel registered data set, where
        # we are grid registered.  In order to make this work
        # we need to adjust the edges of our grid out by half a cell width
        # in each direction.
        txmin = xmin - dx/2.0
        tymax = ymax + dy/2.0

        outshape = (ny, nx)
        transform = Affine.from_gdal(txmin, dx, 0.0, tymax, 0.0, -dy)
        allTouched = not mustContainCenter
        img = features.rasterize(shapes, out_shape=outshape, fill=fillValue,
                                 transform=transform,
                                 all_touched=allTouched,
                                 default_value=burnValue)
        # geodict = GeoDict({'xmin':xmin,'xmax':xmax,'ymin':ymin,'ymax':ymax,
        # 'dx':dx,'dy':dy,'ny':ny,'nx':nx})
        # gd = geodict.asDict()
        # ny,nx = img.shape
        # gd['nx'] = nx
        # gd['ny'] = ny
        # geodict = GeoDict(gd,adjust='bounds')
        return cls(img, geodict)

    def project(self, projection, method='bilinear'):
        """Project Grid2D data into desired projection.

        :param projection:
          Valid proj4 projection string.
        :param method:
          One of the sampling methods described here:
              https://mapbox.github.io/rasterio/topics/resampling.html#resampling-methods
        :raises DataSetException:
          If input projection is not a valid Proj4 string.
          If method is not a valid resampling method found in above URL.
        :returns:
          Re-projected Grid2D object.
        """
        # hack to handle projections that wrap around the 180 meridian.
        crosses = self._geodict.xmax < self._geodict.xmin
        lon_set = 'lon_0' in projection
        old_projection = projection
        if crosses and lon_set:
            # make a new geodict that keeps latitudes the same
            # but centers around lon 0
            tdict = self._geodict.asDict()
            txrange = ((tdict['xmax']+360) - tdict['xmin'])
            tdict['xmin'] = 0 - txrange/2.0
            tdict['xmax'] = 0 + txrange/2.0
            geodict = GeoDict(tdict)

            # make a new projection string centered on lon 0
            projection = _center_projection(projection)
        else:
            geodict = self._geodict

        # check to see if the input projection is valid
        srs = osr.SpatialReference()
        srs.ImportFromProj4(projection)
        if srs.ExportToProj4() == '':
            raise DataSetException('%s is not a valid proj4 string.' %
                                   self._geodict['projection'])

        # check to see if the input resampling method is valid
        int_method = 1  # bi-linear
        try:
            int_method = getattr(Resampling, method)
        except AttributeError:
            raise DataSetException('%s is not a valid resampling '
                                   'method.' % method)

        # get the dimensions of the input data
        nrows, ncols = self._data.shape
        # define the input Affine object
        src_transform = Affine.from_gdal(geodict.xmin -
                                         geodict.dx/2.0,
                                         geodict.dx,
                                         0.0,  # x rotation, not used by us
                                         geodict.ymax +
                                         geodict.dy/2.0,
                                         0.0,  # y rotation, not used by us
                                         # their dy is negative
                                         -1*geodict.dy)

        # set the source and destination projections (have to be CRS
        # dictionaries)
        src_crs = CRS().from_string(geodict.projection).to_dict()
        dst_crs = CRS().from_string(projection).to_dict()

        # determine the boundaries in src coordinates
        if geodict.xmin < geodict.xmax:
            right = geodict.xmax - (geodict.dx/2.0)
        else:
            txmax = geodict.xmax + 360
            right = txmax - (geodict.dx/2.0)
        left = geodict.xmin - (geodict.dx/2.0)
        top = geodict.ymax + (geodict.dy/2.0)
        bottom = geodict.ymin + (geodict.dy/2.0)

        # use this convenience function to determine optimal output
        # transform and dimensions
        dst_transform, width, height = calculate_default_transform(
            src_crs, dst_crs, ncols, nrows, left, bottom, right, top)

        # allocate space for output data (very C-like)
        destination = np.zeros((height, width))

        # if the input has nan values, then tell reproject about that
        # and set the output to that value as well
        src_nan = None
        dst_nan = None
        if np.any(np.isnan(self._data)):
            src_nan = np.nan
            dst_nan = np.nan
        if self._data.dtype in (np.float32, np.float64):
            src_nan = np.nan
            dst_nan = np.nan

        # call the reproject function
        reproject(
            self._data,
            destination,
            src_transform=src_transform,
            src_crs=src_crs,
            dst_transform=dst_transform,
            src_nodata=src_nan,
            dst_nodata=dst_nan,
            dst_crs=projection,
            resampling=int_method)

        # get the pieces of the output transformation
        xmin, dx, xrot, ymax, yrot, mdy = dst_transform.to_gdal()

        # affine dy is negative, so we have to flip it back
        dy = -1*mdy

        # correct for different pixel offsets
        xmin = xmin + (dx/2.0)
        ymax = ymax - (dy/2.0)

        # if we crossed the meridian, we have to set the projection string
        # to reflect where we actually are.
        if crosses and lon_set:
            projection = old_projection

        # Construct a new GeoDict
        gdict = {'xmin': xmin,
                 'xmax': xmin+width*dx,
                 'ymin': ymax-height*dy,
                 'ymax': ymax,
                 'dx': dx,
                 'dy': dy,
                 'nx': width,
                 'ny': height,
                 'projection': projection}
        geodict = GeoDict(gdict, adjust='bounds')

        # Make a new Grid2D object and return it
        newgrid = Grid2D(destination, geodict)
        return newgrid
