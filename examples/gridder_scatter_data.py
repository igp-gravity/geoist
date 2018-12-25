# -*- coding: utf-8 -*-
"""
Created on Mon Dec 17 17:36:09 2018
gridder for scattered data
@author: chens
"""

import pandas as pd
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from geoist.pfm import grdio

#读入数据接口
df1=pd.read_csv("D:\MyWorks\geoist\dataset\gravity-yx.csv")
lon=df1['longitude']
lat=df1['latitude']
grav=df1['elevation']

coordinates = (lon.values, lat.values)

xMin = min(lon.values)
xMax = max(lon.values) 
yMin = min(lat.values)
yMax = max(lat.values)
xNum = 100j
yNum = 80j
nullvalue = 1.701410009187828e+38

grid_x, grid_y = np.mgrid[xMin:xMax:xNum, yMin:yMax:yNum]
#grid_z0 = griddata(coordinates, grav.values, (grid_x, grid_y), method='nearest'
#                   ,fill_value = 0.0 )
#grid_z1 = griddata(coordinates, grav.values, (grid_x, grid_y), method='linear'
#                   ,fill_value = 0.0 )
grid_z2 = griddata(coordinates, grav.values, (grid_x, grid_y), method='cubic'
                   ,fill_value = nullvalue )
#plt.imshow(grid_z0.T)
#plt.imshow(grid_z1.T)
#plt.imshow(grid_z2.T)
grd1=grdio.grddata()
grd1.cols = 100
grd1.rows = 80
grd1.xmin = xMin
grd1.xdim = (xMax-xMin)/(grd1.cols-1)
grd1.ymin = yMin
grd1.ydim = (yMax-yMin)/(grd1.rows-1)

grd1.data = np.ma.masked_equal(grid_z2.T, nullvalue)
grd1.export_surfer(r'D:\MyWorks\geoist\dataset\gravity-yx.grd') 



