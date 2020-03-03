# Name        : dataset.py
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
  DataSet 
"""
# -*- coding: utf-8 -*-
#python 3 compatibility
from __future__ import print_function

#stdlib imports
import abc

#third party imports

#local imports

class DataSetException(Exception):
    """
    Class to represent errors in the DataSet class.
    """
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class DataSetWarning(Warning):
    """
    Class to represent warnings in the DataSet class.
    """
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class DataSet(object):
        #This should be a @classmethod in subclasses
    @abc.abstractmethod
    def load(filename,bounds=None,resample=False,padValue=None):
        """
        Load data into a Grid subclass.  Parameters below are suggested for subclasses.

        :param filename:
           File where data is stored
        :param bounds: 
           Optional tuple of (lonmin,lonmax,latmin,latmax) used to subset data from file.
        :param resample:
           If subsetting data, True indicates that *exact* bounds are desired and data should be resampled to fit.
        :param padValue:
           If asking for data outside bounds of grid, any value not None means fill in those cells with padValue.
           None means don't pad the grid at all.
        :raises NotImplementedError:
           Always for this base class.
        """
        raise NotImplementedError('load method not implemented in base class')

    #TODO: Figure out how format-specific attributes will be handled (ShakeMap, for example)
    #This should be a @classmethod in subclasses
    @abc.abstractmethod
    def save(self,filename): #would we ever want to save a subset of the data?
        """
        Save the data contained in the grid to a format specific file.  Other attributes may be required for
        format specific files.

        :param filename:
           Where file containing data should be written.
        """
        raise NotImplementedError('Save method not implemented in base class')

    @abc.abstractmethod
    def getData(self,getCopy=False):
        """
        Return a reference to or copy of the data inside the Grid

        :param getCopy:
           True indicates that the user wants a copy of the data, not a reference to it.
        :returns:
          A reference to or copy of a numpy array of data.
        """
        raise NotImplementedError('getData method not implemented in base class')

    @abc.abstractmethod
    def setData(self,data):
        """
        Modify the data inside the Grid.

        :param data:
           numpy array of desired data.
        """
        raise NotImplementedError('setData method not implemented in base class')

    @abc.abstractmethod
    def getBounds(self):
        """
        Return the lon/lat range of the data.
        
        :returns:
           Tuple of (lonmin,lonmax,latmin,latmax)
        """
        raise NotImplementedError('getBounds method not implemented in base class')

    @abc.abstractmethod
    def trim(self,geodict,resample=False,method='linear'):
        """
        Trim data to a smaller set of bounds, resampling if requested.  If not resampling,
        data will be trimmed to smallest grid boundary possible.
        
        :param geodict:
           GeoDict object used to specify subset bounds and resolution (if resample is selected)
        :param resample:
           Boolean indicating whether the data should be resampled to *exactly* match input bounds.
        :param method:
           If resampling, method used, one of ('linear','nearest','cubic','quintic')
        """
        raise NotImplementedError('trim method not implemented in base class')

    @abc.abstractmethod
    def getValue(self,lat,lon,method='nearest',default=None): #return nearest neighbor value
        """Return numpy array at given latitude and longitude (using given resampling method).
        
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
        """
        raise NotImplementedError('getValue method not implemented in base class')

    @abc.abstractmethod
    def interpolateToGrid(self,geodict,method='linear'):
        """
        Given a geodict specifying a grid extent and resolution, resample current data set to match.
        
        :param geodict: 
            geodict object from a grid whose extents are inside the extent of this grid.
        :param method: 
            Optional interpolation method - ['linear', 'cubic','quintic','nearest']
        :returns:
          Interpolated grid.
        :raises DataSetException: 
           If the Grid object upon which this function is being called is not completely contained by the grid to which this Grid is being resampled.
        :raises DataSetException: 
           If the resulting interpolated grid shape does not match input geodict.

        This function modifies the internal griddata and geodict object variables.
        """
        raise NotImplementedError('interpolateToGrid method not implemented in base class')
