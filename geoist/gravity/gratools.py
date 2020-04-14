# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 22:28:17 2019

@Author:  Shi CHEN   @ IGP-CEA
         Bei ZHANGE @ IGP-CEA

Co-Author: Jiancang ZHUANG @ ISM-Tokyo

####################################################
            MIT license
        Copyright @ pyBACGS 2019 
         All right reserved 
####################################################

Module which contains the main classes of the pyBACGS
            
CLASS list:
    Outlier 

"""

import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import pygmt as pg
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.io.shapereader import Reader
import shapely.geometry as sgeom
from matplotlib.offsetbox import AnchoredText
from geoist import EXAMPLES_PATH
import geoist.gravity.adjmethods as adj
from geoist import catalog
import geoist.gravity.graobj as ggo


# # Make a plot of the decimated data using Cartopy
# plt.figure(figsize=(7, 6))
# ax = plt.axes(projection=ccrs.Mercator())
# ax.set_title("10' Block Median Bathymetry")
# # Plot the bathymetry as colored circles.
# plt.scatter(lon, lat, c=bathymetry, s=5, transform=ccrs.PlateCarree())
# plt.colorbar().set_label("meters")
# # Use a utility function to setup the tick labels and land feature
# fetch_data.setup_baja_bathymetry_map(ax)
# plt.tight_layout()
# plt.show()

# # Make a plot of the masked data and the data locations.
# crs = ccrs.PlateCarree()
# plt.figure(figsize=(7, 6))
# ax = plt.axes(projection=ccrs.Mercator())
# ax.set_title("Only keep grid points that are close to data")
# ax.plot(data.longitude, data.latitude, ".y", markersize=0.5, transform=crs)
# ax.pcolormesh(*coordinates, dummy_data, transform=crs)
# fetch_data.setup_baja_bathymetry_map(ax, land=None)
# plt.tight_layout()
# plt.show()

class Resviz(object):
    """residual viz
    """
    def __init__(self, res, method):
        """
        """
    pass

class Mapviz(object):
    """
    map viz
    """
    __region = {'cn': [72, 137, 15, 55], 'hb':[106, 125.5, 27, 45.5],
                'hn':[104.5, 123.5, 17, 34], 'nb':[95, 108, 19.5, 44.5],
                'db':[112, 135, 36.5, 54.5], 'xj':[72, 99, 33, 50]}

    def __init__(self, loc = 'cn', engine = 'gmt'):
        """
         region dict = cn, hb, hn, nb, db, xj
        """
        self.region = self.__region[loc]
        self.engine = engine
        
        if engine == 'gmt':
            self.fig = pg.Figure()
            self.prj = "M10i"
        else:
            self.prj = ccrs.PlateCarree()
            self.fig = plt.axes(projection=self.prj)
            
    def plot_base(self, bous = 1, bmap = None, 
                  features = ['land','ocean','lake','river'], 
                  gridline = True, lic = True):
        """
        Plot map base

        Parameters
        ----------
        bous : TYPE, optional
            DESCRIPTION. The default is 1.
        bmap : TYPE, optional
            DESCRIPTION. The default is None.
        features : TYPE, optional
            DESCRIPTION. The default is ['land','ocean','lake','river'].
        gridline : TYPE, optional
            DESCRIPTION. The default is True.
        lic : TYPE, optional
            DESCRIPTION. The default is True.

        Returns
        -------
        None.

        """
        if self.engine == 'gmt':
            self.gmt_plot_base(self.region)
        else:
            self.fig.set_extent([72, 137, 10, 55])
            self.fig.stock_img()
            self.fig.coastlines()
            self.fig.add_feature(cfeature.LAND)
            self.fig.add_feature(cfeature.OCEAN)
            self.fig.add_feature(cfeature.LAKES)
            self.fig.add_feature(cfeature.RIVERS)
            self.fig.gridlines(crs=ccrs.PlateCarree(), draw_labels=True)
            datapath = Path(Path(catalog.__file__).parent,'data')
            fname = Path(datapath, 'bou1_4l.shp')
            f2name = Path(datapath, 'bou2_4l.shp')
            faults = Path(datapath, 'gem_active_faults.shp')
            
            self.fig.add_geometries(Reader(str(faults)).geometries(),
                              ccrs.PlateCarree(),facecolor = 'none',
                              edgecolor='red')
            
            self.fig.add_geometries(Reader(str(f2name)).geometries(),
                              ccrs.PlateCarree(),  facecolor = 'none', 
                              edgecolor='gray', linestyle=':')
            
            self.fig.add_geometries(Reader(str(fname)).geometries(),
                              ccrs.PlateCarree(),  facecolor = 'none', 
                              edgecolor='black')            
                
              
    def plot_points(self):
        pass
        
    def plot_grid(self):
        pass
        
    def gmt_plot_base(self, region, prj = "M10i", frame = True, sls = True, filename = None):
        """
        Using pyGMT plot map     
        Returns
        -------
        None.

        """
        self.fig.basemap(region=region, projection=prj, frame=frame)
        self.fig.coast(shorelines=sls, land="#666666", water="skyblue")


    def gmt_plot_xy(self, data, prj = "M8i", filename = ''):
        """
        Using pyGMT plot xy map     
        Returns
        -------
        None.

        """
        self.fig.plot(
            x=data.longitude,
            y=data.latitude,
            sizes=0.01 * 2 ** data.mag,
            color=data.depth / data.depth.max(),
            cmap="viridis",
            style="cc",
            pen="black",
        )


    def gmt_plot_grd(self, region, grid, prj = "M10i", filename = ''):
        """
        Using pyGMT plot map  
            fig.contour
            fig.grdcontour
            fig.grdimage
        Returns
        -------
        None.

        """
        self.fig.grdimage(grid)


            
    def savefigure(self, filename):
        """
         export gmt figure to disk filename

        Parameters
        ----------
        filename : string

        Returns
        -------
        None.

        """
        if filename == None:
            self.fig.show()
        else:
            self.fig.savefig(filename)
            
class Modanalysis(object):
    pass

class Outlier(object):
    pass

class Gindex(object):
    pass

if __name__ == '__main__':
    
    Mapviz.gmt_plot_base([75, 135, 15, 55], filename = "d:\mapviz1.png")

