# -*- coding: utf-8 -*-
"""
Created on Thu Jan 10 18:34:07 2019
计算WMM2015模型，WMM.cof文件需要放到与py相同目录
@author: chens
"""

import numpy as np
import datetime
from matplotlib.pyplot import show
from geoist.others.scidates import datetime2yeardec
from geoist.pfm.cof import noaamm 

dt = datetime.datetime(2012, 7, 12, 12)
print(datetime2yeardec(dt))
mag = noaamm.noaa(45.5, 105.6, 0.2, datetime2yeardec(dt), mod='wmm')
print("F:",mag.total.item()) #F
print("D:",mag.decl.item())  #D
print("I:",mag.incl.item())  #I

lon, lat = np.meshgrid(np.arange(-180, 180+10, 10), np.arange(-90, 90+10, 10))
mag = noaamm.noaa(lat, lon, 0, 2015)
noaamm.plotwmm(mag)
