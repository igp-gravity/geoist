# -*- coding: utf-8 -*-
"""
Created on Thu Jan 10 18:34:07 2019
计算IGRF12模型，需要IGRF*.cof文件
计算WMM2015模型，WMM**.cof文件需要放到与py相同目录,cof文件在geoist.pfm.cof目录下
@author: chens
"""

import numpy as np
import datetime
from geoist.others.scidates import datetime2yeardec
from geoist.pfm.cof import noaamm 
import geoist.pfm.igrf as igrf

dt = datetime.datetime(2017, 7, 12, 12)
print('Date：', datetime2yeardec(dt))
#IGRF 12 model
m1=igrf.IGRF()
mag_igrf=m1.pnt(45.5, 105.6, 0.2, datetime2yeardec(dt))
print('IGRF Model(F I D):\n',mag_igrf)

#WMM 2015 model
print('WMM2015 Model(F I D):')
mag = noaamm.noaa(45.5, 105.6, -1.0, datetime2yeardec(dt), mod='wmm')
print("F:",mag.total.item()) #F
print("I:",mag.incl.item())  #I
print("D:",mag.decl.item())  #D


lon, lat = np.meshgrid(np.arange(-180, 180+10, 10), np.arange(-90, 90+10, 10))
mag = noaamm.noaa(lat, lon, 0, 2015)
noaamm.plotwmm(mag)




