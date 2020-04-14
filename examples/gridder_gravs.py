# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 13:40:07 2020

@author: chens
"""

import matplotlib.pyplot as plt
from geoist.others.gdal import GDALGrid
from geoist.others.utils import map2DGrid,Grid2xyz

import warnings
warnings.simplefilter("ignore")

filename1 = 'D:\\research-work\\Topex\\WGM2012\\WGM2012_ETOPO1_ponc_2min.grd'
filename2 = 'D:\\research-work\\Topex\\WGM2012\\WGM2012_Freeair_ponc_2min.grd'
filename3 = 'D:\\research-work\\Topex\\WGM2012\\WGM2012_Bouguer_ponc_2min.grd'
filename4 = 'D:\\research-work\\Topex\\WGM2012\\WGM2012_Isostatic_ponc_2min.grd'

# 通过GDAL库读取网格数据
gd1, ff = GDALGrid.getFileGeoDict(filename1)
gd2, ff = GDALGrid.getFileGeoDict(filename2)
gd3, ff = GDALGrid.getFileGeoDict(filename3)

print(gd1)
gd1.xmin = 70.0
gd1.xmax = 105.0
gd1.ymin = 15.0
gd1.ymax = 45.0
gd1.nx = 800
gd1.ny = 600

# 读取指定区域
topo = GDALGrid.load(filename1, samplegeodict = gd1, resample = True)
fag = GDALGrid.load(filename2, samplegeodict = gd1, resample = True)
bug = GDALGrid.load(filename3, samplegeodict = gd1, resample = True)
iso = GDALGrid.load(filename4, samplegeodict = gd1, resample = True)

# 取数据
topoc = topo.cut(80,90, 30,35, align = True)
tcx, tcy, tcz = Grid2xyz(topoc) # 位场变换要求x,y,z三列数据格式接口 pftrans

p_jw = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
p_merc = "+proj=merc +lon_0=100 +ellps=WGS84 +datum=WGS84 +no_defs +units=km"
p_lcc = "+proj=lcc +lon_0=100 +lat_0=35 +lat_1=45 +ellps=WGS84 +datum=WGS84 +no_defs +units=km"

# 投影
topo1 = topo.project(p_merc)
topo2 = topo.project(p_lcc)

fig, (ax0, ax1) = plt.subplots(nrows=1, ncols=2,figsize=(12,6))
ax0.imshow(topo1.getData())
ax1.imshow(topo2.getData())


topo0 = topo2.project(p_jw)
plt.figure()
plt.imshow(topo0.getData())

# 画图

fig, axes = plt.subplots(nrows=2,ncols=2,figsize=(12,12))
fig.tight_layout()
map2DGrid(axes[0,0],topo,'Topography',10,10, isLeft=True)
map2DGrid(axes[0,1],fag,'Free-air gravity',10,10)
map2DGrid(axes[1,0],bug,'Bouguer gravity',10,10, isLeft=True)
map2DGrid(axes[1,1],iso,'Isostatic gravity',10,10)
