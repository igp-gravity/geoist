# -*- coding: utf-8 -*-
"""
 Name        : gmt.py
 Created on  : 2018/09/11 17:00
 Author      : Steve Chen<chenshi@cea-igp.ac.cn>
 Affiliation : Institute of Geophysics, CEA.
 Version     : 0.1.0
 Copyright   : Copyright (C) 2018-2020 GEOIST Development Team. All Rights Reserved.
 License     : Distributed under the MIT License. See LICENSE.txt for more info.
 Github      : https://igp-gravity.github.io/
 Description : Application for GMT io
"""

#python 3 compatibility
from __future__ import print_function

#stdlib imports
import struct
import os.path
import sys

#third party imports
import numpy as np
from scipy.io import netcdf
from .grid2d import Grid2D
from .gridbase import Grid
from .dataset import DataSetException
from .geodict import GeoDict
import h5py



'''Grid2D subclass for reading,writing, and manipulating GMT format grids.

Usage:

::

     gmtgrid = GMTGrid.load(gmtfilename)
     gmtgrid.getGeoDict()

This class supports reading and writing of all three GMT formats: NetCDF, HDF, and the GMT "native" format.

'''

NETCDF_TYPES = {'B':np.uint8,
                'b':np.int8,
                'h':np.int16,
                'i':np.int32,
                'f':np.float32,
                'd':np.float64}

INVERSE_NETCDF_TYPES = {'uint8':'B',
                        'int8':'b',
                        'int16':'h',
                        'int32':'i',
                        'float32':'f',
                        'float64':'d'}

def subsetArray(data,data_range,fgeodict):
    iulx1 = data_range['iulx1']
    iuly1 = data_range['iuly1']
    ilrx1 = data_range['ilrx1']
    ilry1 = data_range['ilry1']
    data1 = data[iuly1:ilry1,iulx1:ilrx1]
    ymax1,xmin1 = fgeodict.getLatLon(iuly1,iulx1)
    ymin1,xmax1 = fgeodict.getLatLon(ilry1-1,ilrx1-1)
    if 'iulx2' in data_range:
        iulx2 = data_range['iulx2']
        iuly2 = data_range['iuly2']
        ilrx2 = data_range['ilrx2']
        ilry2 = data_range['ilry2']
        data2 = data[iuly2:ilry2,iulx2:ilrx2]
        data = np.hstack((data1,data2)).copy()
        ymax2,xmin2 = fgeodict.getLatLon(iuly2,iulx2)
        ymin2,xmax2 = fgeodict.getLatLon(ilry2-1,ilrx2-1)
        ny,nx = data.shape
        ymin2,xmax2 = fgeodict.getLatLon(ilry2-1,ilrx2-1)
        xmax = xmax2
    else:
        data = data1.copy()
        xmax = xmax1
    ny,nx = data.shape
    geodict = GeoDict({'xmin':xmin1,
                       'xmax':xmax,
                       'ymin':ymin1,
                       'ymax':ymax1,
                       'dx':fgeodict.dx,
                       'dy':fgeodict.dy,
                       'nx':nx,
                       'ny':ny})
    return (data,geodict)

def sub2ind(shape,subtpl):
    """
    Convert 2D subscripts into 1D index.
    @param shape: Tuple indicating size of 2D array.
    @param subtpl: Tuple of (possibly) numpy arrays of row,col values.
    @return: 1D array of indices.
    """
    if len(shape) != 2 or len(shape) != len(subtpl):
        raise IndexError("Input size and subscripts must have length 2 and be equal in length")
    
    row,col = subtpl
    ny,nx = shape
    ind = nx*row + col
    return ind

def indexArray(array,shp,i1,i2,j1,j2):
    if not isinstance(j1,int) and len(j1):
        j1 = j1[0]
    if not isinstance(j2,int) and len(j2):
        j2 = j2[0]
    if not isinstance(i1,int) and len(i1):
        i1 = i1[0]
    if not isinstance(i2,int) and len(i2):
        i2 = i2[0]
    if len(array.shape) == 1:
        ny = i2-i1
        nx = j2-j1
        if hasattr(array,'dtype'):
            data = np.zeros((ny,nx),dtype=array.dtype)
        else:
            typecode = array.typecode()
            dtype = NETCDF_TYPES[typecode]
            data = np.zeros((ny,nx),dtype=dtype)
        rowidx = np.arange(i1,i2)
        i = 0
        for row in rowidx:
            idx1 = sub2ind(shp,(row,j1))
            idx2 = sub2ind(shp,(row,j2))
            data[i,:] = array[idx1:idx2]
            i += 1
    else:
        ny,nx = array.shape
        i1r = ny-i1
        i2r = ny-i2
        data = array[i2r:i1r,j1:j2].copy()
    return data

def createSampleXRange(M,N,filename,bounds=None,dx=None,dy=None):
    if dx is None:
        dx = 1.0
    if dy is None:
        dy = 1.0
    if bounds is None:
        xmin = 0.5
        xmax = xmin + (N-1)*dx
        ymin = 0.5
        ymax = ymin + (M-1)*dy
    else:
        xmin,xmax,ymin,ymax = bounds
    data = np.arange(0,M*N).reshape(M,N).astype(np.int32)
    cdf = netcdf.netcdf_file(filename,'w')
    cdf.createDimension('side',2)
    cdf.createDimension('xysize',M*N)
    dim = cdf.createVariable('dimension','i',('side',))
    dim[:] = np.array([N,M])
    spacing = cdf.createVariable('spacing','i',('side',))
    spacing[:] = np.array([dx,dy])
    zrange = cdf.createVariable('z_range',INVERSE_NETCDF_TYPES[str(data.dtype)],('side',))
    zrange[:] = np.array([data.min(),data.max()])
    x_range = cdf.createVariable('x_range','d',('side',))
    x_range[:] = np.array([xmin,xmax])
    y_range = cdf.createVariable('y_range','d',('side',))
    y_range[:] = np.array([ymin,ymax])
    z = cdf.createVariable('z',INVERSE_NETCDF_TYPES[str(data.dtype)],('xysize',))
    z[:] = data.flatten()
    cdf.close()
    return data

def createSampleGrid(M,N):
    """Used for internal testing - create an NxN grid with lower left corner at 0.5,0.5, dx/dy = 1.0.
    :param M:
       Number of rows in output grid
    :param N:
       Number of columns in output grid
    :returns:
       GMTGrid object where data values are an NxN array of values from 0 to N-squared minus 1, and geodict
       lower left corner is at 0.5/0.5 and cell dimensions are 1.0.
    """
    data = np.arange(0,M*N).reshape(M,N)
    data = data.astype(np.int32) #arange gives int64 by default, not supported by netcdf3
    xvar = np.arange(0.5,0.5+N,1.0)
    yvar = np.arange(0.5,0.5+M,1.0)
    geodict = {'ny':N,
               'nx':N,
               'xmin':0.5,
               'xmax':xvar[-1],
               'ymin':0.5,
               'ymax':yvar[-1],
               'dx':1.0,
               'dy':1.0}
    gmtgrid = GMTGrid(data,geodict)
    return gmtgrid

class GMTGrid(Grid2D):
    '''Grid2D subclass for reading,writing, and manipulating GMT format grids.

    Usage:

    ::

     gmtgrid = GMTGrid.load(gmtfilename)
     gmtgrid.getGeoDict()

     This class supports reading and writing of all three GMT formats: NetCDF, HDF, and the GMT "native" format.

     '''
    def __init__(self,data,geodict):
        """Construct a GMTGrid object.
        :param data:
           2D numpy data array (must match geodict spec)
        :param geodict:
           Dictionary specifying the spatial extent,resolution and shape of the data.
        :returns:
           A GMTGrid object.
        :raises DataSetException:
          When data and geodict dimensions do not match. 
        """
        m,n = data.shape
        if m != geodict.ny or n != geodict.nx:
            raise DataSetException('Input geodict does not match shape of input data.')
        self._data = data
        self._geodict = geodict

    @classmethod
    def getFileType(cls,grdfile):
        """Get the GMT sub-format (netcdf, hdf, or GMT binary).
        :param grdfile:
          File name of suspected GMT grid file.
        :returns:
          One of 'netcdf' (NetCDF version 3), 'hdf' (NetCDF version 4), 'native' (so-called GMT native format), or 'unknown'.
        """
        isThree = True
        if sys.version_info.major == 2:
            isThree = False
        f = open(grdfile,'rb')
        #check to see if it's HDF or CDF
        ftype = 'unknown'
        try:
            f.seek(1,0)
            if isThree:
                hdfsig = f.read(3).decode('utf-8')
            else:
                hdfsig = ''.join(struct.unpack('ccc',f.read(3)))
            
            if hdfsig == 'HDF':
                ftype = 'hdf'
            else:
                f.seek(0,0)
                if isThree:
                    cdfsig = f.read(3).decode('utf-8')
                else:
                    cdfsig = ''.join(struct.unpack('ccc',f.read(3)))
                if cdfsig == 'CDF':
                    ftype = 'netcdf'
                else:
                    f.seek(0,0)
                    nx = struct.unpack('I',f.read(4))[0]
                    f.seek(8,0)
                    offset = struct.unpack('I',f.read(4))[0]
                    if (offset == 0 or offset == 1) and nx > 0:
                        ftype = 'native'
        except:
            pass            
        f.close()
        return ftype

    @classmethod
    def getFileGeoDict(cls,filename):
        """Get the spatial extent, resolution, and shape of grid inside NetCDF file.
        :param filename:
           File name of NetCDF file.
        :returns:
          GeoDict specifying spatial extent, resolution, and shape of grid inside NetCDF file.
        :raises DataSetException:
          When the file is not detectable as one of the three flavors of GMT grids.
        """
        ftype = cls.getFileType(filename)
        if ftype == 'native':
            geodict,xvar,yvar,fmt,zscale,zoffset = cls.getNativeHeader(filename)
        elif ftype == 'netcdf':
            geodict,xvar,yvar = cls.getNetCDFHeader(filename)
        elif ftype == 'hdf':
            geodict,xvar,yvar = cls.getHDFHeader(filename)
        else:
            raise DataSetException('Unknown file type for file "%s".' % filename)
        newgeodict,first_column_duplicated = cls.checkFirstColumnDuplicated(geodict)
        return (newgeodict,first_column_duplicated)
    
    @classmethod
    def getNativeHeader(cls,fname,fmt=None):
        """Get the header information from a GMT native grid file.
        :param fname:
           File name of GMT native grid
        :param fmt:
          One of:
            - 'h' for 16 bit signed integer
            - 'i' for 32 bit signed integer
            - 'f' for 32 bit float
            - 'd' for 64 bit float
        :returns:
           - GeoDict specifying spatial extent, resolution, and shape of grid inside NetCDF file.
           - xvar array specifying X coordinates of data columns
           - yvar array specifying Y coordinates of data rows
           - fmt If input fmt is None, this will be a best guess as to the data format in the file.
           - zscale Data multiplier
           - zoffset Value to be added to data
        """
        f = open(fname,'rb')
        #Given that we can't automatically distinguish between 32 bit ints and 32 bit floats, we'll use
        #this value as a cutoff for "detecting" when a value read from disk as a float has an "unreasonably" 
        #high exponent value.  This is NOT guaranteed to work - use the fmt keyword if you want to be sure.
        MAX_FLOAT_EXP = 30
        fsize = os.path.getsize(fname)
        datalen = fsize-892
        f.seek(0,0)
        geodict = {}
        geodict['nx'] = struct.unpack('I',f.read(4))[0]
        geodict['ny'] = struct.unpack('I',f.read(4))[0]
        offset = struct.unpack('I',f.read(4))[0]
        geodict['xmin'] = struct.unpack('d',f.read(8))[0]
        geodict['xmax'] = struct.unpack('d',f.read(8))[0]

        #GMT grids can have longitudes from 0-360, we're correcting for that here
        if geodict['xmax'] > 180:
            geodict['xmax'] = geodict['xmax'] - 360
        if geodict['xmax'] > 180:
            geodict['xmin'] = geodict['xmin'] - 360
            
        geodict['ymin'] = struct.unpack('d',f.read(8))[0]
        geodict['ymax'] = struct.unpack('d',f.read(8))[0]
        zmin = struct.unpack('d',f.read(8))[0]
        zmax = struct.unpack('d',f.read(8))[0]
        geodict['dx'] = struct.unpack('d',f.read(8))[0]
        geodict['dy'] = struct.unpack('d',f.read(8))[0]
        zscale = struct.unpack('d',f.read(8))[0]
        zoffset = struct.unpack('d',f.read(8))[0]
        xunits = f.read(80).strip()
        yunits = f.read(80).strip()
        zunits = f.read(80).strip()
        title = f.read(80).strip()
        command = f.read(320).strip()
        remark = f.read(160).strip()
        npixels = geodict['ny']*geodict['nx']
        lenshort = npixels*2
        lenfloat = npixels*4
        lendouble = npixels*8
        if fmt is None:
            if datalen == lenshort:
                fmt = 'h'
            if datalen == lendouble:
                fmt = 'd'
            if datalen == lenfloat: #let's try to guess whether this is float or int
                fpos = f.tell()
                #read 1 byte, check to see if it's nan or 0 - if it is, then we definitely have a float
                dbytes = struct.unpack('f',f.read(4))[0]
                while dbytes == 0.0:
                    dbytes = struct.unpack('f',f.read(4))[0]
                f.seek(fpos) #go back to where we were
                if np.isnan(dbytes):
                    fmt = 'f'
                elif int(np.abs(np.log10(dbytes))) > MAX_FLOAT_EXP: #does this have a crazy large exponent?
                    fmt = 'i'
                else:
                    fmt = 'f'
        f.close()

        #We are going to represent all grids internally as grid-line registered
        #The difference between pixel and gridline-registered grids is depicted well here:
        #http://gmt.soest.hawaii.edu/doc/5.1.0/GMT_Docs.html#grid-registration-the-r-option
        if offset == 1:
            geodict['xmin'] += geodict['dx']/2.0
            geodict['xmax'] -= geodict['dx']/2.0
            geodict['ymin'] += geodict['dy']/2.0
            geodict['ymax'] -= geodict['dy']/2.0
        xvar,dx2 = np.linspace(geodict['xmin'],geodict['xmax'],num=geodict['nx'],retstep=True)
        yvar,dy2 = np.linspace(geodict['ymin'],geodict['ymax'],num=geodict['ny'],retstep=True)
        gd = GeoDict(geodict)
        return (gd,xvar,yvar,fmt,zscale,zoffset)

    @classmethod
    def readFile(cls,filename,data_range):
        """
        Read any GMT formatted file.
        :param filename:
          Input GMT formatted grid file.
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
        :returns:
           data is a 2D numpy array of all data found inside bounds, and 
          geodict gives the geo-referencing information for the data.
        """
        ftype = cls.getFileType(filename)
        if ftype == 'netcdf':
            data,geodict = cls.readNetCDF(filename,data_range)
        elif ftype == 'hdf':
            data,geodict = cls.readHDF(filename,data_range)
        elif ftype == 'native':
            data,geodict = cls.readGMTNative(filename,data_range)
        return (data,geodict)
    
    @classmethod
    def readNetCDF(cls,filename,data_range):
        """Read a NetCDF formatted GMT file.
        
        :param filename:
          Input GMT formatted grid file.
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
        :returns:
           data is a 2D numpy array of all data found inside bounds, and 
          geodict gives the geo-referencing information for the data.
        """
        fgeodict,first_column_duplicated = cls.getFileGeoDict(filename)
        dx = fgeodict.dx
        dy = fgeodict.dy
        gny = fgeodict.ny
        gnx = fgeodict.nx
        iulx1 = data_range['iulx1']
        iuly1 = data_range['iuly1']
        ilrx1 = data_range['ilrx1']
        ilry1 = data_range['ilry1']
        cdf = netcdf.netcdf_file(filename)
        if 'z' in cdf.variables:
            zvar = cdf.variables['z']
        else:
            zvar = cdf.variables['elevation']
        isScanLine = len(zvar.shape) == 1
        if isScanLine:
            data = indexArray(zvar,(gny,gnx),0,gny,0,gnx)
        else:
            data = np.flipud(indexArray(zvar,(gny,gnx),0,gny,0,gnx))
        if first_column_duplicated:
            data = data[:,0:-1]

        newdata,geodict = subsetArray(data,data_range,fgeodict)
        del data
        del zvar
        cdf.close()
        
        return (newdata,geodict)

    @classmethod
    def readGMTNative(cls,fname,data_range,fmt=None):
        """Read the data and geo-referencing information from a GMT native grid file, subsetting if requested.
        
        http://gmt.soest.hawaii.edu/doc/5.1.2/GMT_Docs.html#native-binary-grid-files

        :param fname:
          File name of GMT native grid
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
        :param fmt: 
           Data width, one of: 
             - 'i' (16 bit signed integer)
             - 'l' (32 bit signed integer)
             - 'f' (32 bit float)
             - 'd' (64 bit float)

          Strictly speaking, this is only necessary when the data file
          is 32 bit float or 32 bit integer, as there is no *sure* way
          to tell from the header or data which data type is contained
          in the file.  If fmt is None, then the code will try to
          guess as best it can from the data whether it is integer or
          floating point data. Caveat emptor!

          :returns:
          Tuple of data (2D numpy array of data, possibly subsetted from file) and geodict (see above).
        """
        HDRLEN = 892
        fgeodict,xvar,yvar,fmt,zscale,zoffset = cls.getNativeHeader(fname,fmt)
        t,first_column_duplicated = cls.getFileGeoDict(fname)
        sfmt = '%i%s' % (fgeodict.nx*fgeodict.ny,fmt)
        dwidths = {'h':2,'i':4,'f':4,'d':8}
        dwidth = dwidths[fmt]
        f = open(fname,'rb')
        f.seek(HDRLEN)
        #TODO - do something clever with memory mapping some day so that
        #we don't always have to read the whole file just to subset it.
        #Currently not an issue.
        dbytes = f.read(fgeodict.nx*fgeodict.ny*dwidth)
        f.close()

        data = np.flipud(np.array(struct.unpack(sfmt,dbytes)))
        data.shape = (fgeodict.ny,fgeodict.nx)
        if zscale != 1.0 or zoffset != 0.0:
            data = (data * zscale) + zoffset
        if first_column_duplicated:
            data = data[:,0:-1]

        data,geodict = subsetArray(data,data_range,fgeodict)

        return (data,geodict)        

    @classmethod
    def getNetCDFHeader(cls,filename):
        """Get the header information from a GMT NetCDF3 file.
        :param fname:
           File name of GMT NetCDF3 grid
        :returns:
           - GeoDict specifying spatial extent, resolution, and shape of grid inside NetCDF file.
           - xvar array specifying X coordinates of data columns
           - xvar array specifying Y coordinates of data rows
        """
        cdf = netcdf.netcdf_file(filename)
        geodict = {}
        xvarname = None
        registration = 'gridline'
        if hasattr(cdf,'node_offset') and getattr(cdf,'node_offset') == 1:
            registration = 'pixel'
        if 'x' in cdf.variables.keys():
            xvarname = 'x'
            yvarname = 'y'
        elif 'lon' in cdf.variables.keys():
            xvarname = 'lon'
            yvarname = 'lat'
        if xvarname is not None:
            xvar = cdf.variables[xvarname].data.copy()

            #GMT grids can have longitudes from 0-360, we're correcting for that here
            #however, sometimes GMT (or at least netcdf) allows an array for longitude centers
            #where the max value exceeds +180.0.  Here we're also trapping for that situation, 
            #and hopefully not catastrophically breaking a bunch of other stuff.
            if (xvar > 180.01).any():
                xvar = xvar - 360
                
            yvar = cdf.variables[yvarname].data.copy()
            geodict['nx'] = len(xvar)
            geodict['ny'] = len(yvar)
            geodict['xmin'] = xvar.min()
            geodict['xmax'] = xvar.max()
            geodict['ymin'] = yvar.min()
            geodict['ymax'] = yvar.max()
            newx = np.linspace(geodict['xmin'],geodict['xmax'],num=geodict['nx'])
            newy = np.linspace(geodict['ymin'],geodict['ymax'],num=geodict['ny'])
            geodict['dx'] = newx[1]-newx[0]
            geodict['dy'] = newy[1]-newy[0]
        elif 'x_range' in cdf.variables.keys():
            geodict['xmin'] = cdf.variables['x_range'].data[0]
            geodict['xmax'] = cdf.variables['x_range'].data[1]
            geodict['ymin'] = cdf.variables['y_range'].data[0]
            geodict['ymax'] = cdf.variables['y_range'].data[1]
            geodict['nx'],geodict['ny'] = cdf.variables['dimension'].data
            #geodict['dx'],geodict['dy'] = cdf.variables['spacing'].data
            xvar = np.linspace(geodict['xmin'],geodict['xmax'],num=geodict['nx'])
            yvar = np.linspace(geodict['ymin'],geodict['ymax'],num=geodict['ny'])
            geodict['dx'] = xvar[1] - xvar[0]
            geodict['dy'] = yvar[1] - yvar[0]
        else:
            raise DataSetException('No support for CDF data file with variables: %s' % str(cdf.variables.keys()))

        #We are going to represent all grids internally as grid-line registered
        #The difference between pixel and gridline-registered grids is depicted well here:
        #http://gmt.soest.hawaii.edu/doc/5.1.0/GMT_Docs.html#grid-registration-the-r-option
        if registration == 'pixel':
            geodict['xmin'] += geodict['dx']/2.0
            geodict['xmax'] -= geodict['dx']/2.0
            geodict['ymin'] += geodict['dy']/2.0
            geodict['ymax'] -= geodict['dy']/2.0
        #because dx/dy are not explicitly defined in netcdf headers, here we'll assume
        #that those values are adjustable, and we'll preserve the shape and extent.
        gd = GeoDict(geodict,adjust='res')
        return (gd,xvar,yvar)

    @classmethod
    def _subsetRegions(self,zvar,sampledict,fgeodict,xvar,yvar,firstColumnDuplicated):
        """Internal method used to do subsampling of data for all three GMT formats.
        :param zvar:
          A numpy array-like thing (CDF/HDF variable, or actual numpy array)
        :param sampledict:
          GeoDict indicating the bounds where data should be sampled.
        :param fgeodict:
          Geo dictionary with the file information.
        :param xvar:
          Numpy array specifying X coordinates of data columns
        :param yvar:
          Numpy array specifying Y coordinates of data rows
        :param firstColumnDuplicated:
          Boolean - is this a file where the last column of data is the same as the first (for grids that span entire globe).
        :returns:
          Tuple of (data,geodict) (subsetted data and geodict describing that data).
        """
        isScanLine = len(zvar.shape) == 1
        txmin,txmax,tymin,tymax = (sampledict.xmin,sampledict.xmax,sampledict.ymin,sampledict.ymax)
        #we're not doing anything fancy with the data here, just cutting out what we need
        if fgeodict.xmin > fgeodict.xmax:
            fxmax = fgeodict.xmax + 360
        else:
            fxmax = fgeodict.xmax
        xmin = max(fgeodict.xmin,txmin)
        xmax = min(fxmax,txmax)
        ymin = max(fgeodict.ymin,tymin)
        ymax = min(fgeodict.ymax,tymax)
        #these are the bounds of the whole file
        gxmin = fgeodict.xmin
        gxmax = fgeodict.xmax
        gymin = fgeodict.ymin
        gymax = fgeodict.ymax
        dx = fgeodict.dx
        dy = fgeodict.dy
        gny = fgeodict.ny
        gnx = fgeodict.nx
        geodict = fgeodict.copy().asDict()
        if xmin == gxmin and xmax == gxmax and ymin == gymin and ymax == gymax:
            #read in the whole file
            #data = np.flipud(zvar[:].copy())
            if isScanLine:
                data = indexArray(zvar,(gny,gnx),0,gny,0,gnx)
            else:
                data = np.flipud(indexArray(zvar,(gny,gnx),0,gny,0,gnx))
            if firstColumnDuplicated:
                data = data[:,0:-1]
                geodict['xmax'] -= geodict['dx']
        else:
            #what are the nearest grid coordinates to our desired bounds?
            #for example, if the grid starts at xmin = 0.5 and xmax = 6.5 with dx=1.0, 
            #and the user wants txmin = 2.0 and txmax of 4.0, the pixel coordinates
            #(erring on the side of including more data), would be at txmin = 1.5 and 
            #txmax = 4.5.
            txmin2 = gxmin + dx*np.floor((txmin - gxmin)/dx)
            txmax2 = gxmin + dx*np.ceil((txmax - gxmin)/dx)
            tymin2 = gymin + dy*np.floor((tymin - gymin)/dy)
            tymax2 = gymin + dy*np.ceil((tymax - gymin)/dy)
            if xmin > xmax:
                #cut user's request into two regions - one from the minimum to the
                #meridian, then another from the meridian to the maximum.
                iuly1,iulx1 = fgeodict.getRowCol(tymax2,txmin2)
                ilry1,ilrx1 = fgeodict.getRowCol(tymin2,fgeodict.xmax)
                #get section from the 180 meridian to xmax
                iuly2,iulx2 = fgeodict.getRowCol(tymax2,fgeodict.xmin)
                ilry2,ilrx2 = fgeodict.getRowCol(tymin2,txmax2)

                #remove a column if the first column is duplicated in the data set
                if firstColumnDuplicated:
                    ilrx1 -= 1

                tny = (ilry1 - iuly1)+1
                tnx = (ilrx1 - iulx1)+1 + (ilrx2 - iulx2)+1
                
                #in Python 3 the "long" type was integrated into the int type.
                if sys.version_info.major == 2:
                    outcols1 = long(ilrx1-iulx1)
                    outcols2 = long(ilrx2-iulx2)
                    outcols = long(outcols1+outcols2)
                    outrows = long(ilry1-iuly1)
                else:
                    outcols1 = int(ilrx1-iulx1)
                    outcols2 = int(ilrx2-iulx2)
                    outcols = int(outcols1+outcols2)
                    outrows = int(ilry1-iuly1)
                section1 = indexArray(zvar,(gny,gnx),iuly1,ilry1+1,iulx1,ilrx1+1)
                #section1 = zvar[iuly1:ilry1,iulx1:ilrx1].copy()
                #section2 = zvar[iuly2:ilry2,iulx2:ilrx2].copy()
                section2 = indexArray(zvar,(gny,gnx),iuly2,ilry2+1,iulx2,ilrx2+1)
                if isScanLine:
                    data = np.concatenate((section1,section2),axis=1)
                else:
                    data = np.flipud(np.concatenate((section1,section2),axis=1))
                outrows,outcols = data.shape
                newymax,newxmin = fgeodict.getLatLon(iuly1,iulx1)
                newymin,newxmax = fgeodict.getLatLon(ilry2,ilrx2)
                geodict['xmin'] = newxmin
                geodict['xmax'] = newxmax
                geodict['ymin'] = newymin
                geodict['ymax'] = newymax
                geodict['ny'],geodict['nx'] = data.shape
            else:
                iuly,iulx = fgeodict.getRowCol(tymax2,txmin2)
                ilry,ilrx = fgeodict.getRowCol(tymin2,txmax2)
                tny = (ilry - iuly)+1
                tnx = (ilrx - iulx)+1
                    
                if isScanLine:
                    data = indexArray(zvar,(gny,gnx),iuly,ilry+1,iulx,ilrx+1)
                else:
                    data = np.flipud(indexArray(zvar,(gny,gnx),iuly,ilry+1,iulx,ilrx+1))

                newymax,newxmin = fgeodict.getLatLon(iuly,iulx)
                newymin,newxmax = fgeodict.getLatLon(ilry,ilrx)
                geodict['xmin'] = newxmin
                geodict['xmax'] = newxmax
                geodict['ymin'] = newymin
                geodict['ymax'] = newymax
                geodict['ny'],geodict['nx'] = data.shape
                
        gd = GeoDict(geodict)
        return (data,gd)
    
    @classmethod
    def getHDFHeader(cls,hdffile):
        """Get the header information from a GMT NetCDF4 (HDF) file.
        :param fname:
           File name of GMT NetCDF4 grid
        :returns:
          GeoDict specifying spatial extent, resolution, and shape of grid inside NetCDF file.
        """
        geodict = {}
        f = h5py.File(hdffile,'r')
        registration = 'gridline'
        if f.get('node_offset') is not None and f.attrs['node_offset'][0] == 1:
            registration = 'pixel'
        if 'x' in f.keys():
            xvarname = 'x'
            yvarname = 'y'
        elif 'lon' in f.keys():
            xvarname = 'lon'
            yvarname = 'lat'
        if xvarname is not None:
            xvar = f[xvarname][:]
            yvar = f[yvarname][:]

            #GMT grids can have longitudes from 0-360, we're correcting for that here
            #however, sometimes GMT (or at least hdf) allows an array for longitude centers
            #where the max value exceeds +180.0.  Here we're also trapping for that situation, 
            #and hopefully not catastrophically breaking a bunch of other stuff.
            if (xvar > 180.01).any():
                xvar = xvar - 360
            
            geodict['ny'] = len(yvar)
            geodict['nx'] = len(xvar)
            geodict['xmin'] = xvar[0]
            geodict['xmax'] = xvar[-1]
            geodict['ymin'] = yvar[0]
            geodict['ymax'] = yvar[-1]
            newx = np.linspace(geodict['xmin'],geodict['xmax'],num=geodict['nx'])
            newy = np.linspace(geodict['ymin'],geodict['ymax'],num=geodict['ny'])
            geodict['dx'] = newx[1]-newx[0]
            geodict['dy'] = newy[1]-newy[0]
        else:
            geodict['xmin'] = f['x_range'][0]
            geodict['xmax'] = f['x_range'][1]
            geodict['ymin'] = f['y_range'][0]
            geodict['ymax'] = f['y_range'][1]
            geodict['nx'],geodict['ny'] = (f['dimension'][0],f['dimension'][1])
            xvar = np.linspace(geodict['xmin'],geodict['xmax'],num=nx)
            yvar = np.linspace(geodict['ymin'],geodict['ymax'],num=ny)
            geodict['dx'] = xvar[1] - xvar[0]
            geodict['dy'] = yvar[1] - yvar[0]

        #We are going to represent all grids internally as grid-line registered
        #The difference between pixel and gridline-registered grids is depicted well here:
        #http://gmt.soest.hawaii.edu/doc/5.1.0/GMT_Docs.html#grid-registration-the-r-option
        if registration == 'pixel':
            geodict['xmin'] += geodict['dx']/2.0
            geodict['xmax'] -= geodict['dx']/2.0
            geodict['ymin'] += geodict['dy']/2.0
            geodict['ymax'] -= geodict['dy']/2.0
        f.close()
        #because dx/dy are not explicitly defined in hdf headers, here we'll assume
        #that those values are adjustable, and we'll preserve the shape and extent.
        gd = GeoDict(geodict,adjust='res')
        return (gd,xvar,yvar)

    @classmethod
    def readHDF(cls,hdffile,data_range):
        """Read a NetCDF formatted GMT file.
        
        :param filename:
          Input GMT formatted grid file.
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
        :returns:
           data is a 2D numpy array of all data found inside bounds, and 
          geodict gives the geo-referencing information for the data.
        """
        #fgeodict,xvar,yvar = cls.getHDFHeader(hdffile)
        fgeodict,first_column_duplicated = cls.getFileGeoDict(hdffile)
        f = h5py.File(hdffile,'r')
        zvar = f['z']
        data = np.flipud(zvar[:])
        if first_column_duplicated:
            data = data[:,0:-1]
        data,geodict = subsetArray(data,data_range,fgeodict)
        f.close()
        return (data,geodict)
    
    def save(self,filename,format='netcdf'):
        """Save a GMTGrid object to a file.
        :param filename:
          Name of desired output file.
        :param format:
          One of 'netcdf','hdf' or 'native'.
        :raises DataSetException:
          When format not one of ('netcdf,'hdf','native')
        """
        if format not in ['netcdf','hdf','native']:
            raise DataSetException('Only NetCDF3, HDF (NetCDF4), and GMT native output are supported.')
        if format == 'netcdf':
            f = netcdf.NetCDFFile(filename,'w')
            m,n = self._data.shape
            dx = f.createDimension('x',n)
            dy = f.createDimension('y',m)
            x = f.createVariable('x',np.float64,('x'))
            y = f.createVariable('y',np.float64,('y'))
            x[:] = np.linspace(self._geodict.xmin,self._geodict.xmax,self._geodict.nx)
            y[:] = np.linspace(self._geodict.ymin,self._geodict.ymax,self._geodict.ny)
            z = f.createVariable('z',self._data.dtype,('y','x'))
            z[:] = np.flipud(self._data)
            z.actual_range = np.array((np.nanmin(self._data),np.nanmax(self._data)))
            f.close()
        elif format == 'hdf':
            #Create the file and the top-level attributes GMT wants to see
            f = h5py.File(filename,'w')
            f.attrs['Conventions'] = 'COARDS, CF-1.5'
            f.attrs['title'] = 'filename'
            f.attrs['history'] = 'Created with python GMTGrid.save(%s,format="hdf")' % filename
            f.attrs['GMT_version'] = 'NA'

            #Create the x array and the attributes of that GMT wants to see
            xvar = np.linspace(self._geodict.xmin,self._geodict.xmax,self._geodict.nx)
            x = f.create_dataset('x',data=xvar,shape=xvar.shape,dtype=str(xvar.dtype))
            x.attrs['CLASS'] = 'DIMENSION_SCALE'
            x.attrs['NAME'] = 'x'
            x.attrs['_Netcdf4Dimid'] = 0 #no idea what this is
            x.attrs['long_name'] = 'x'
            x.attrs['actual_range'] = np.array((xvar[0],xvar[-1]))

            #Create the x array and the attributes of that GMT wants to see
            yvar = np.linspace(self._geodict.ymin,self._geodict.ymax,self._geodict.ny)
            y = f.create_dataset('y',data=yvar,shape=yvar.shape,dtype=str(yvar.dtype))
            y.attrs['CLASS'] = 'DIMENSION_SCALE'
            y.attrs['NAME'] = 'y'
            y.attrs['_Netcdf4Dimid'] = 1 #no idea what this is
            y.attrs['long_name'] = 'y'
            y.attrs['actual_range'] = np.array((yvar[0],yvar[-1]))
            
            #create the z data set
            z = f.create_dataset('z',data=np.flipud(self._data),shape=self._data.shape,dtype=str(self._data.dtype))
            z.attrs['long_name'] = 'z'
            #zvar.attrs['_FillValue'] = array([ nan], dtype=float32)
            z.attrs['actual_range'] = np.array((np.nanmin(self._data),np.nanmax(self._data)))
            
            #close the hdf file
            f.close()
        elif format == 'native':
            f = open(filename,'wb')
            f.write(struct.pack('I',self._geodict.nx))
            f.write(struct.pack('I',self._geodict.ny))
            f.write(struct.pack('I',0)) #gridline registration
            f.write(struct.pack('d',self._geodict.xmin))
            f.write(struct.pack('d',self._geodict.xmax))
            f.write(struct.pack('d',self._geodict.ymin))
            f.write(struct.pack('d',self._geodict.ymax))
            f.write(struct.pack('d',self._data.min()))
            f.write(struct.pack('d',self._data.max()))
            f.write(struct.pack('d',self._geodict.dx))
            f.write(struct.pack('d',self._geodict.dy))
            f.write(struct.pack('d',1.0)) #scale factor to multiply data by
            f.write(struct.pack('d',0.0)) #offfset to add to data
            f.write(struct.pack('80s',b'X units (probably degrees)'))
            f.write(struct.pack('80s',b'Y units (probably degrees)'))
            f.write(struct.pack('80s',b'Z units unknown'))
            f.write(struct.pack('80s',b'')) #title
            f.write(struct.pack('320s',b'Created with GMTGrid() class, a product of the NEIC.')) #command
            f.write(struct.pack('160s',b'')) #remark
            if self._data.dtype not in [np.int16,np.int32,np.float32,np.float64]:
                msg = 'Data type of "%s" is not supported by the GMT native format.'
                raise DataSetException(msg % str(self._data.dtype))
            fpos1 = f.tell()
            #the left-right flip is necessary because of the way tofile() works
            newdata = np.fliplr(np.flipud(self._data[:])) 
            newdata.tofile(f)
            fpos2 = f.tell()
            bytesout = fpos2 - fpos1
            f.close()
    
            
    @classmethod
    def load(cls,filename,samplegeodict=None,resample=False,method='linear',doPadding=False,padValue=np.nan):
        """Create a GMTGrid object from a (possibly subsetted, resampled, or padded) GMT grid file.
        :param gmtfilename:
          Name of input file.
        :param samplegeodict:
          GeoDict used to specify subset bounds and resolution (if resample is selected)
        :param resample:
          Boolean used to indicate whether grid should be resampled from the file based on samplegeodict.
        :param method:
          If resample=True, resampling method to use ('nearest','linear','cubic','quintic')
        :param doPadding:
          Boolean used to indicate whether, if samplegeodict is outside bounds of grid, to pad values around the edges.
        :param padValue:
          Value to fill in around the edges if doPadding=True.
        :returns:
          GMTgrid instance (possibly subsetted, padded, or resampled)
        :raises DataSetException:
          * When sample bounds are outside (or too close to outside) the bounds of the grid and doPadding=False.
          * When the input file type is not recognized.
        """
        #get the geodict describing the source file, plus a boolean telling us if the last column
        #is a duplicate of the first column
        filegeodict,first_column_duplicated = cls.getFileGeoDict(filename)

        #if the sampling geodict is identical to the file geodict, then turn resampling off
        if samplegeodict is not None and samplegeodict == filegeodict:
            resample = False

        #If the sample grid is aligned with the host grid, then resampling won't accomplish anything 
        # if samplegeodict is not None and filegeodict.isAligned(samplegeodict):
        #     resample = False
        
        #buffer out the sample geodict (if resampling) enough to allow interpolation.
        if samplegeodict is not None:
            sampledict = cls.bufferBounds(samplegeodict,filegeodict,resample=resample,doPadding=doPadding) #parent static method
        else:
            sampledict = filegeodict
            
        #Ensure that the two grids at least 1) intersect and 2) are aligned if resampling is True.
        try:
            cls.verifyBounds(filegeodict,sampledict,resample=resample) #parent static method, may raise an exception
        #if they're not, and we have padding on, then give them a grid with all pad values.
        except DataSetException as dse:
            if doPadding:
                if not filegeodict.contains(sampledict):
                    data = np.ones((sampledict.ny,sampledict.nx),dtype=np.float32)*padValue
                    return cls(data=data,geodict=sampledict)
                
        sampledict = filegeodict.getIntersection(sampledict)
        bounds = (sampledict.xmin,sampledict.xmax,sampledict.ymin,sampledict.ymax)
        
        data_range = cls.getDataRange(filegeodict,sampledict,
                                      first_column_duplicated=first_column_duplicated)
        data,geodict = cls.readFile(filename,data_range)
        pad_dict = cls.getPadding(filegeodict,samplegeodict,doPadding=doPadding) #parent static method
        data,geodict = cls.padGrid(data,geodict,pad_dict)
        grid = cls(data=data,geodict=geodict)
        if resample:
            grid = grid.interpolateToGrid(samplegeodict,method=method)
            
        if np.any(np.isinf(grid._data)):
            grid._data[np.isinf(grid._data)] = padValue
        return grid
            
        
class BinCDFArray(object):
    def __init__(self,array,ny,nx):
        self.array = array
        self.ny = ny
        self.nx = nx

    def __getitem__(self,*args):
        """Allows slicing of CDF data array in the same way as a numpy array."""
        if len(args) == 1 and isinstance(args[0][0],int):
            #user has passed in a tuple of row,col - they only want one value
            row = args[0][0]
            col = args[0][1]
            ny = self.ny
            nx = self.nx
            if row < 0 or row > ny-1:
                raise Exception("Row index out of bounds")
            if col < 0 or col > nx-1:
                raise Exception("Row index out of bounds")
            idx = nx * row + col
            offset = 0
            return self.array[idx]

        if len(args) == 1 and isinstance(args[0][0],slice): #they want a non-scalar subset of the data
            ny = self.ny
            nx = self.nx
            #calculate offset to first data element
            key1 = args[0][0]
            key2 = args[0][1]
            rowstart = key1.start
            rowend = key1.stop
            rowstep = key1.step
            colstart = key2.start
            colend = key2.stop
            colstep = key2.step
            
            if rowstep is None:
                rowstep = 1
            if colstep is None:
                colstep = 1

            #error checking
            if rowstart < 0 or rowstart > ny-1:
                raise Exception("Row index out of bounds")
            if rowend < 0 or rowend > ny:
                raise Exception("Row index out of bounds")
            if colstart < 0 or colstart > nx-1:
                raise Exception("Col index out of bounds")
            if colend < 0 or colend > nx:
                raise Exception("Col index out of bounds")

            colcount = (colend-colstart)
            rowcount = (rowend-rowstart)
            outrows = np.ceil(rowcount/rowstep)
            outcols = np.ceil(colcount/colstep)
            data = np.zeros([outrows,outcols],dtype=self.dtype)
            outrow = 0
            for row in range(int(rowstart),int(rowend),int(rowstep)):
                #just go to the beginning of the row, we're going to read in the whole line
                idx = nx*row 
                offset = self.dwidth*idx #beginning of row
                line = self.array[idx:idx+nx]
                data[outrow,:] = line[colstart:colend:colstep]
                outrow = outrow+1
                
        else:
            raise Exception("Unsupported __getitem__ input %s" % (str(key)))
        return(data)

        
        
