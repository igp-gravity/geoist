# -*- coding: utf-8 -*-
"""
Created on Sun Dec  2 11:36:59 2018
使用giplt画等值线图
@author: chens
"""
from geoist import gridder
import matplotlib.pyplot as plt
from geoist.vis import giplt

area = [-20, 20, -50, 50]
x, y = gridder.scatter(area, n=100)
data = x**2 + y**2
plt.figure()
plt.axis('scaled')
giplt.contourf(y, x, data, shape=(50, 50), levels=30, interp=True)
plt.colorbar(orientation='horizontal')
plt.plot(y, x, '.k')
plt.xlabel('y (East-West)')
plt.ylabel('x (North-South)')
plt.show()
