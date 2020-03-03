# Name        : gdal.py
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
"""
 Description : functions for GDAL drivers
"""
# -*- coding: utf-8 -*-
#stdlib imports
import os.path
import sys
from collections import OrderedDict
import warnings

#third party imports
import rasterio
from affine import Affine
import numpy as np
from .grid2d import Grid2D
from .dataset import DataSetException,DataSetWarning
from .geodict import GeoDict

def get_affine(src):
    aff = None
    # See https://github.com/mapbox/rasterio/issues/86
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        aff = src.transform
    if isinstance(aff,list):
        aff = Affine.from_gdal(*aff)
    return aff

class GDALGrid(Grid2D):
    def __init__(self,data,geodict):
        """Construct a GMTGrid object.
        :param data:
           2D numpy data array (must match geodict spec)
        :param geodict:
           GeoDict Object specifying the spatial extent,resolution and shape of the data.
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
    def load(cls,filename,samplegeodict=None,resample=False,method='linear',doPadding=False,padValue=np.nan):
        """
        This method should do the following:
        1) If resampling, buffer bounds outwards.
        2) Translate geographic bounds to pixel bounds (Grid2D static method)
        3) Create a Grid2D subclass instance
        4) Call readFile() with pixel bounds
        5) Pad as requested.
        6) Resample as requested.
        7) Return Grid2D subclass.
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

        if samplegeodict is not None and not filegeodict.intersects(samplegeodict):
             if not doPadding:
                 raise DataSetException('Your sampling grid is not contained by the file.  To load anyway, use doPadding=True.')
        
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

    @classmethod
    def getFileGeoDict(cls,filename):
        """Get the spatial extent, resolution, and shape of grid inside ESRI grid file.
        :param filename:
           File name of ESRI grid file.
        :returns:
           - GeoDict object specifying spatial extent, resolution, and shape of grid inside ESRI grid file.
        :raises DataSetException:
          When the file contains a grid with more than one band.
          When the file geodict is internally inconsistent.
        """
        geodict = {}

        with rasterio.open(filename) as src:
            aff = get_affine(src)
            if aff is None:
                fmt = 'Could not find .transform attribute from GDAL dataset.'
                raise AttributeError(fmt)
            
            geodict['dx'] = aff.a
            geodict['dy'] = -1*aff.e
            geodict['xmin'] = aff.xoff + geodict['dx']/2.0
            geodict['ymax'] = aff.yoff - geodict['dy']/2.0
            

            shp = src.shape
            if len(shp) > 2:
                raise DataSetException('Cannot support grids with more than one band')
            geodict['ny'] = src.height
            geodict['nx'] = src.width
            geodict['xmax'] = geodict['xmin'] + (geodict['nx']-1)*geodict['dx']
            if geodict['xmax'] == geodict['xmin']:
                pass

            geodict['ymin'] = geodict['ymax'] - (geodict['ny']-1)*geodict['dy']

            gd = GeoDict(geodict)

        newgeodict,first_column_duplicated = cls.checkFirstColumnDuplicated(gd)
        return (newgeodict,first_column_duplicated)
    
    @classmethod
    def _subsetRegions(self,src,sampledict,fgeodict,firstColumnDuplicated):
        """Internal method used to do subsampling of data for all three GMT formats.
        :param zvar:
          A numpy array-like thing (CDF/HDF variable, or actual numpy array)
        :param sampledict:
          GeoDict object with bounds and row/col information.
        :param fgeodict:
          GeoDict object with the file information.
        :param firstColumnDuplicated:
          Boolean - is this a file where the last column of data is the same as the first (for grids that span entire globe).
        :returns:
          Tuple of (data,geodict) (subsetted data and geodict describing that data).
        """
        txmin,txmax,tymin,tymax = (sampledict.xmin,sampledict.xmax,sampledict.ymin,sampledict.ymax)
        trows,tcols = (sampledict.ny,sampledict.nx)
        if fgeodict.xmin > fgeodict.xmax:
            fxmax = fgeodict.xmax + 360
        else:
            fxmax = fgeodict.xmax
        #we're not doing anything fancy with the data here, just cutting out what we need
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
        geodict = None
        if xmin == gxmin and xmax == gxmax and ymin == gymin and ymax == gymax:
            #just read the whole file
            tfdict = fgeodict.asDict()
            data = src.read()
            data = np.squeeze(data)
            if firstColumnDuplicated:
                data = data[:,0:-1]
                tfdict['xmax'] -= geodict.dx
                tfdict['nx'] -= 1
            geodict = GeoDict(tfdict)
        else:
            #what are the nearest grid coordinates to our desired bounds?
            #for example, if the grid starts at xmin = 0.5 and xmax = 6.5 with dx=1.0, 
            #and the user wants txmin = 2.0 and txmax of 4.0, the pixel coordinates
            #(erring on the side of including more data), would be at txmin = 1.5 and 
            #txmax = 4.5.
            if not fgeodict.isAligned(sampledict):
                txmin2 = gxmin + dx*np.floor((txmin - gxmin)/dx)
                txmax2 = gxmin + dx*np.ceil((txmax - gxmin)/dx)
                tymin2 = gymin + dy*np.floor((tymin - gymin)/dy)
                tymax2 = gymin + dy*np.ceil((tymax - gymin)/dy)
            else:
                txmin2 = xmin
                txmax2 = xmax
                tymin2 = ymin
                tymax2 = ymax
            if txmin2 > txmax2:
                #cut user's request into two regions - one from the minimum to the
                #meridian, then another from the meridian to the maximum.
                #new create sections algorithm
                #get section from the xmin to the 180 meridian
                iuly1,iulx1 = fgeodict.getRowCol(tymax2,txmin2)
                ilry1,ilrx1 = fgeodict.getRowCol(tymin2,fgeodict.xmax)
                #get section from the 180 meridian to xmax
                iuly2,iulx2 = fgeodict.getRowCol(tymax2,fgeodict.xmin)
                ilry2,ilrx2 = fgeodict.getRowCol(tymin2,txmax2)

                if firstColumnDuplicated:
                    ilrx1 -= 1

                tny = (ilry1 - iuly1)+1
                tnx = (ilrx1 - iulx1)+1 + (ilrx2 - iulx2)+1
                
                #(region1,region2) = self._createSections((xmin,xmax,ymin,ymax),fgeodict,firstColumnDuplicated)
                #(iulx1,iuly1,ilrx1,ilry1) = region1
                #(iulx2,iuly2,ilrx2,ilry2) = region2
                window1 = ((iuly1,ilry1+1),(iulx1,ilrx1+1))
                window2 = ((iuly2,ilry2+1),(iulx2,ilrx2+1))
                section1 = src.read(1,window=window1)
                section2 = src.read(1,window=window2)
                data = np.hstack((section1,section2))
                tfdict = {}
                newymax,newxmin = fgeodict.getLatLon(iuly1,iulx1)
                newymin,newxmax = fgeodict.getLatLon(ilry2,ilrx2)
                tfdict['xmin'] = newxmin
                tfdict['xmax'] = newxmax
                tfdict['ymin'] = newymin
                tfdict['ymax'] = newymax
                tfdict['dx'] = dx
                tfdict['dy'] = dy
                tfdict['ny'],tfdict['nx'] = data.shape
                geodict = GeoDict(tfdict)
            else:
                iuly,iulx = fgeodict.getRowCol(tymax2,txmin2)
                ilry,ilrx = fgeodict.getRowCol(tymin2,txmax2)
                tny = (ilry - iuly)+1
                tnx = (ilrx - iulx)+1

                window = ((iuly,ilry+1),(iulx,ilrx+1))
                tfdict = {}
                newymax,newxmin = fgeodict.getLatLon(iuly,iulx)
                newymin,newxmax = fgeodict.getLatLon(ilry,ilrx)
                tfdict['xmin'] = newxmin
                tfdict['xmax'] = newxmax
                tfdict['ymin'] = newymin
                tfdict['ymax'] = newymax
                tfdict['dx'] = dx
                tfdict['dy'] = dy
                #window = ((iymin,iymax+1),(ixmin,ixmax+1))
                data = src.read(1,window=window)
                data = np.squeeze(data)
                tfdict['ny'],tfdict['nx'] = data.shape
                geodict = GeoDict(tfdict)
            
        
        return (data,geodict)

    @classmethod
    def readFile(cls,filename,data_range):
        """
        Read an ESRI flt/bip/bil/bsq formatted file using rasterIO (GDAL Python wrapper).
        :param filename:
          Input ESRI formatted grid file.
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
          A tuple of (data,geodict) where data is a 2D numpy array of all data found inside bounds, and 
          geodict gives the geo-referencing information for the data.
        """
        iulx1 = data_range['iulx1']
        ilrx1 = data_range['ilrx1']
        iuly1 = data_range['iuly1']
        ilry1 = data_range['ilry1']
        if 'iulx2' in data_range:
            iulx2 = data_range['iulx2']
            ilrx2 = data_range['ilrx2']
            iuly2 = data_range['iuly2'] 
            ilry2 = data_range['ilry2']
        else:
            iulx2 = None
            ilrx2 = None
            iuly2 = None
            ilry2 = None
        data = None

        with rasterio.open(filename) as src:
            window1 = ((iuly1,ilry1),
                       (iulx1,ilrx1))
            section1 = src.read(1,window=window1)
            if 'iulx2' in data_range:
                window2 = ((iuly2,ilry2),
                           (iulx2,ilrx2))
                section2 = src.read(1,window=window2)
                data = np.hstack((section1,section2))
            else:
                data = section1

        #Put NaN's back in where nodata value was
        nodata = src.get_nodatavals()[0]
        if nodata is not None and data.dtype in [np.float32,np.float64]: #NaNs only valid for floating point data
            if (data==nodata).any():
                data[data == nodata] = np.nan

        ny,nx = data.shape
        filegeodict,first_column_duplicated = cls.getFileGeoDict(filename)
        ymax1,xmin1 = filegeodict.getLatLon(iuly1,iulx1)
        ymin1,xmax1 = filegeodict.getLatLon(ilry1-1,ilrx1-1)
        xmin = xmin1
        ymax = ymax1
        ymin = ymin1
        xmax = xmax1
        if iulx2 is not None:
            ymin2,xmax2 = filegeodict.getLatLon(ilry2-1,ilrx2-1)
            xmax = xmax2

        dx = filegeodict.dx
        dy = filegeodict.dy
        geodict = GeoDict({'xmin':xmin,
                           'xmax':xmax,
                           'ymin':ymin,
                           'ymax':ymax,
                           'nx':nx,
                           'ny':ny,
                           'dx':dx,
                           'dy':dy},adjust='res')
        
            
        return (data,geodict)
                

    def _getHeader(self):
        hdr = {}
        if sys.byteorder == 'little':
            hdr['BYTEORDER'] = 'LSBFIRST'
        else:
            hdr['BYTEORDER'] = 'MSBFIRST'
        hdr['LAYOUT'] = 'BIL'
        hdr['NROWS'],hdr['NCOLS'] = self._data.shape
        hdr['NBANDS'] = 1
        if self._data.dtype == np.uint8:
            hdr['NBITS'] = 8
            hdr['PIXELTYPE'] = 'UNSIGNEDINT'
        elif self._data.dtype == np.int8:
            hdr['NBITS'] = 8
            hdr['PIXELTYPE'] = 'SIGNEDINT'
        elif self._data.dtype == np.uint16:
            hdr['NBITS'] = 16
            hdr['PIXELTYPE'] = 'UNSIGNEDINT'
        elif self._data.dtype == np.int16:
            hdr['NBITS'] = 16
            hdr['PIXELTYPE'] = 'SIGNEDINT'
        elif self._data.dtype == np.uint32:
            hdr['NBITS'] = 32
            hdr['PIXELTYPE'] = 'UNSIGNEDINT'
        elif self._data.dtype == np.int32:
            hdr['NBITS'] = 32
            hdr['PIXELTYPE'] = 'SIGNEDINT'
        elif self._data.dtype == np.float32:
            hdr['NBITS'] = 32
            hdr['PIXELTYPE'] = 'FLOAT'
        elif self._data.dtype == np.float64:
            hdr['NBITS'] = 32
            hdr['PIXELTYPE'] = 'FLOAT'
        else:
            raise DataSetException('Data type "%s" not supported.' % str(self._data.dtype))
        hdr['BANDROWBYTES'] = hdr['NCOLS']*(hdr['NBITS']/8)
        hdr['TOTALROWBYTES'] = hdr['NCOLS']*(hdr['NBITS']/8)
        hdr['ULXMAP'] = self._geodict.xmin
        hdr['ULYMAP'] = self._geodict.ymax
        hdr['XDIM'] = self._geodict.dx
        hdr['YDIM'] = self._geodict.dy
        #try to have a nice readable NODATA value in the header file
        zmin = np.nanmin(self._data)
        zmax = np.nanmax(self._data)
        if self._data.dtype in [np.int8,np.int16,np.int32]:
            nodata = np.array([-1*int('9'*i) for i in range(3,20)])
            if zmin > nodata[-1]:
                NODATA = nodata[np.where(nodata < zmin)[0][0]]
            else: #otherwise just pick an arbitrary value smaller than our smallest
                NODATA = zmin - 1
        else:
            nodata = np.array([int('9'*i) for i in range(3,20)])
            if zmin < nodata[-1]:
                NODATA = nodata[np.where(nodata > zmin)[0][0]]
            else: #otherwise just pick an arbitrary value smaller than our smallest
                NODATA = zmax + 1
        hdr['NODATA'] = NODATA
        keys = ['BYTEORDER','LAYOUT','NROWS','NCOLS','NBANDS','NBITS','BANDROWBYTES','TOTALROWBYTES','PIXELTYPE',
                'ULXMAP','ULYMAP','XDIM','YDIM','NODATA']
        hdr2 = OrderedDict()
        for key in keys:
            hdr2[key] = hdr[key]
        return hdr2
    
    def save(self,filename,format='EHdr'):
        """
        Save the data contained in this grid to a float or integer ESRI grid file.  Described here:
        http://webhelp.esri.com/arcgisdesktop/9.3/index.cfm?TopicName=BIL,_BIP,_and_BSQ_raster_files
        http://resources.esri.com/help/9.3/arcgisdesktop/com/gp_toolref/conversion_tools/float_to_raster_conversion_.htm.

        :param filename:
          String representing file to which data should be saved.
        :param format:
          Currently this code only supports the GDAL format 'EHdr' (see formats above.)  As rasterIO write support is expanded, this code should add functionality accordingly.
        :raises DataSetException:
          When format is not 'EHdr'.        
        """
        supported = ['EHdr']
        if format not in supported:
            raise DataSetException('Only "%s" file formats supported for saving' % str(supported))
        hdr = self._getHeader()
        data = self._data #create a reference to the data - this may be overridden by a downcasted version for doubles
        if self._data.dtype == np.float32:
            data = self._data.astype(np.float32) #so we can find/reset nan values without screwing up original data
            data[np.isnan(data)] = hdr['NODATA']
        elif self._data.dtype == np.float64:
            data = self._data.astype(np.float32)
            data[np.isnan(data)] = hdr['NODATA']
            warnings.warn(DataSetWarning('Down-casting double precision floating point to single precision'))

        data.tofile(filename)
        #write out the header file
        basefile,ext = os.path.splitext(filename)
        hdrfile = basefile+'.hdr'
        f = open(hdrfile,'wt')
        for (key,value) in hdr.items():
            value = hdr[key]
            f.write('%s  %s\n' % (key,str(value)))
        f.close()
            
    



