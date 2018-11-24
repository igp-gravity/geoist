# -*- coding: utf-8 -*-
"""
 Name        : gridcontainer.py
 Created on  : 2018/09/11 17:00
 Author      : Steve Chen<chenshi@cea-igp.ac.cn>
 Affiliation : Institute of Geophysics, CEA.
 Version     : 0.1.0
 Copyright   : Copyright (C) 2018-2020 GEOIST Development Team. All Rights Reserved.
 License     : Distributed under the MIT License. See LICENSE.txt for more info.
 Github      : https://igp-gravity.github.io/
 Description : Application for grid processing.
"""

# stdlib imports
from datetime import datetime
import collections
import time
import io
import copy

# third party imports
import h5py
import numpy as np
from impactutils.io.container import HDFContainer

# local imports
from .grid2d import Grid2D
from .geodict import GeoDict

GROUPS = {'grid': 'grids'}


class GridHDFContainer(HDFContainer):
    def setGrid(self, name, grid, metadata=None, compression=True):
        """Store a Grid2D object as a dataset.

        Args:
          name (str): Name of the Grid2D object to be stored.
          grid (Grid2D): Grid2D object to be stored.
          metadata (dict): Simple dictionary (values of strings and numbers).
          compression (bool): Boolean indicating whether dataset should be compressed
                              using the gzip algorithm.

        Returns:
          HDF Dataset containing grid and metadata.
        """
        if compression:
            compression = 'gzip'
        else:
            compression = None

        grid_name = '%s' % name
        if GROUPS['grid'] not in self._hdfobj:
            grid_group = self._hdfobj.create_group(GROUPS['grid'])
        else:
            grid_group = self._hdfobj[GROUPS['grid']]

        array_metadata = grid.getGeoDict().asDict()
        data = grid.getData()
        if metadata is not None:
            array_metadata.update(metadata)
        dset = grid_group.create_dataset(
            grid_name, data=data, compression=compression)
        for key, value in array_metadata.items():
            dset.attrs[key] = value
        return dset

    def getGrid(self, name):
        """
        Retrieve a Grid2D object and any associated metadata from the container.

        Args:
            name (str):
                The name of the Grid2D object stored in the container.

        Returns:
            (tuple) Grid2D object, and a dictionary of metadata.
        """
        grid_name = '%s' % name
        grid_group = self._hdfobj[GROUPS['grid']]
        if grid_name not in grid_group:
            raise LookupError('Grid %s not in %s'
                              % (name, self.getFileName()))

        dset = grid_group[grid_name]
        data = dset[()]
        array_metadata, meta_metadata = _split_dset_attrs(dset)
        geodict = GeoDict(array_metadata)
        grid = Grid2D(data, geodict)
        return grid, meta_metadata

    def getGrids(self):
        """
        Return list of names of Grid2D objects stored in container.

        Returns:
          (list) List of names of Grid2D objects stored in container.
        """
        if GROUPS['grid'] not in self._hdfobj:
            return []
        grids = list(self._hdfobj[GROUPS['grid']].keys())
        return grids

    def dropGrid(self, name):
        """
        Delete Grid2D object from container.

        Args:
          name (str):
                The name of the Grid2D object to be deleted.

        """
        mgrid = '%s' % name
        grid_group = self._hdfobj[GROUPS['grid']]
        if mgrid not in grid_group:
            raise LookupError('Grid %s not in %s'
                              % (name, self._hdfobj.filename))
        del grid_group[mgrid]


def _split_dset_attrs(dset):
    geodict = {}
    metadata = {}
    grid_keys = ['xmin', 'xmax', 'ymin', 'ymax', 'nx', 'ny', 'dx', 'dy']
    for key, value in dset.attrs.items():
        if key in grid_keys:
            geodict[key] = value
        else:
            metadata[key] = value
    return (geodict, metadata)
