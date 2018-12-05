# -*- coding: utf-8 -*-
"""
Created on Sun Dec  2 11:36:59 2018
使用giplt画投影地图
@author: chens
"""
# 3rd packages
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap 
# local packages
from geoist.vis import giplt
from geoist import gridder

area = [-20, 20, -50, 50]
x, y = gridder.scatter(area, n=100)
data = x**2 + y**2

plt.figure()
bm = giplt.basemap(area, projection='robin')
giplt.draw_countries(bm)
giplt.draw_coastlines(bm)
giplt.contourf(y, x, data, shape=(50, 50), levels=30,interp=True, basemap=bm)
plt.colorbar(orientation='horizontal')
plt.show()