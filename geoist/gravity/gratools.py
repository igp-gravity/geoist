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

import geoist.gravity.adjmethods as adj
import geoist.gravity.graobj as ggo
import pygmt as pg
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
    __region = {'cn': [75, 135, 15, 55], 'hb':[75, 135, 15, 55],
                'hn':[75, 135, 15, 55],'nb':[75, 135, 15, 55],
                'db':[75, 135, 15, 55],'xj':[75, 135, 15, 55]}

    def __init__(self, loc = 'cn', engine = 'gmt'):
        """
        region dict = 
        """
        self.region = self.__region[loc]
        self.engine = engine

    @classmethod
    def gmt_plot_base(cls, region, prj = "M10i", frame = True, sls = True, filename = ''):
        """
        Using pyGMT plot map     
        Returns
        -------
        None.

        """
        fig = pg.Figure()
        fig.basemap(region=region, projection=prj, frame=frame)
        fig.coast(shorelines=sls, land="#666666", water="skyblue")
        if (len(filename)>0):
            fig.savefig(filename)
        else:
            fig.show()

    @classmethod
    def gmt_plot_xy(cls, region, data, prj = "M8i", frame = True, sls = True, filename = ''):
        """
        Using pyGMT plot xy map     
        Returns
        -------
        None.

        """
        fig = pg.Figure()
        fig.basemap(region=region, projection=prj, frame=frame)
        fig.coast(shorelines=sls, land="black", water="skyblue")
        fig.plot(
            x=data.longitude,
            y=data.latitude,
            sizes=0.02 * 2 ** data.magnitude,
            color=data.depth_km / data.depth_km.max(),
            cmap="viridis",
            style="cc",
            pen="black",
        )
        if (len(filename)>0):
            fig.savefig(filename)
        else:
            fig.show()

    @classmethod
    def gmt_plot_grd(cls, region, prj = "M10i", frame = True, sls = True, filename = ''):
        """
        Using pyGMT plot map  
            fig.contour
            fig.grdcontour
            fig.grdimage
        Returns
        -------
        None.

        """
        fig = pg.Figure()
        fig.basemap(region=region, projection=prj, frame=frame)
        fig.coast(shorelines=sls)
        if (len(filename)>0):
            fig.savefig(filename)
        else:
            fig.show()

class Modanalysis(object):
    pass

class Outlier(object):
    pass

class Gindex(object):
    pass

if __name__ == '__main__':
    
    Mapviz.gmt_plot_base([75, 135, 15, 55], filename = "d:\mapviz1.png")

