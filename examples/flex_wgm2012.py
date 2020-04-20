# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 13:40:07 2020

@author: chens
"""

import matplotlib.pyplot as plt
from geoist.others.gdal import GDALGrid
from geoist.others.utils import map2DGrid
from geoist.others.grid2d import Grid2D

import warnings
warnings.simplefilter("ignore")

filename1 = 'D:\\research-work\\Topex\\WGM2012\\WGM2012_ETOPO1_ponc_2min.grd'
filename2 = 'D:\\research-work\\Topex\\WGM2012\\WGM2012_Freeair_ponc_2min.grd'
filename3 = 'D:\\research-work\\Topex\\WGM2012\\WGM2012_Bouguer_ponc_2min.grd'

# 通过GDAL库读取网格数据
gd1, ff = GDALGrid.getFileGeoDict(filename1)
gd2, ff = GDALGrid.getFileGeoDict(filename2)
gd3, ff = GDALGrid.getFileGeoDict(filename3)

print(gd1)
gd1.xmin = 70.0
gd1.xmax = 105.0
gd1.ymin = 15.0
gd1.ymax = 45.0
gd1.nx = 400
gd1.ny = 300


# 读取指定区域
topo = GDALGrid.load(filename1, samplegeodict = gd1, resample = True)
fag = GDALGrid.load(filename2, samplegeodict = gd1, resample = True)
bug = GDALGrid.load(filename3, samplegeodict = gd1, resample = True)

# 取数据
p_jw = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
#p_lcc = "+proj=lcc +lon_0=90 +lat_0=35 +lat_1=45 +ellps=WGS84 +datum=WGS84 +no_defs +units=km"
p_lcc = "+proj=lcc +lon_0=90 +lat_0=35 +lat_1=45 +ellps=WGS84 +datum=WGS84 +no_defs +units=m"

# 投影
topo_m = topo.project(p_lcc)
bug_m = bug.project(p_lcc)
#topoc = topo_km.cut(-1000,400, 1000,2000, align = True)

#import geoist.flex as gflex
from geoist.flex import TopoGrid, BougGrid, Project

topoflex = TopoGrid(topo_m.getData()[::-1], topo_m._geodict.dx/1000, topo_m._geodict.dy/1000)
topoflex.filter_water_depth()
bugflex = BougGrid(bug_m.getData()[::-1], bug_m._geodict.dx/1000, bug_m._geodict.dy/1000)

contours = topoflex.make_contours(0.)
mask = (topoflex.data < -500.)

topoflex.plot(mask=mask, contours=contours, cmap='gist_earth', vmin=-1000, vmax=7000)
bugflex.plot(mask=mask, contours=contours, cmap='seismic', vmin=-400, vmax=400)

project = Project(grids=[topoflex, bugflex])
project.init()
project.wlet_admit_coh() #need time
# project.__dict__.keys()
project.plot_admit_coh(kindex=7, contours=contours, mask=mask)
# Take random cell value within grid and set as attribute
project.cell = (250, 100)
# Plot admittance and coherence functions
project.plot_functions()
# # EET
#project.inverse = 'L2'
project.inverse = 'bayes'
print(project.inverse)
project.mask = mask
project.estimate_grid(30, atype='joint')
project.plot_results(mean_Te=True, mask=True, contours=contours, cmap='Spectral')
project.plot_results(std_Te=True, mask=True, contours=contours, cmap='magma_r')


eetdict = topo_m.getGeoDict().copy()
eetny, eetnx = project.mean_Te_grid.shape
eetdict.nx = eetnx
eetdict.ny = eetny
eetgrid = Grid2D(project.mean_Te_grid[::-1],eetdict)
eet_jw = eetgrid.project(p_jw)
# transform to jw 
topo_jw = topo_m.project(p_jw)

fig, (ax0, ax1) = plt.subplots(nrows=1, ncols=2,figsize=(12,6))
ax0.imshow(topo_jw.getData())
ax1.imshow(topo.getData())

fig, axes = plt.subplots(nrows=2,ncols=2,figsize=(16,12))
fig.tight_layout()
map2DGrid(axes[0,0],topo,'Topography',5,5, isLeft=True, xinc = 5.0)
map2DGrid(axes[0,1],topo_jw,'Topography',5,5, isLeft=True, xinc = 5.0)
map2DGrid(axes[1,0],topo,'Topography',5,5, isLeft=True, prj = 'merc')
map2DGrid(axes[1,1],topo_jw,'Topography',5,5, prj = 'merc') #, gmap = 2, cmap ='seismic')

fig, axes = plt.subplots(nrows=1,ncols=1,figsize=(12,12))
map2DGrid(axes,eet_jw,'Elastic thickness',5,5, isLeft=True, xinc = 5.0
          , bous = 1, cmap = 'Spectral' )
# map2DGrid(axes[1],topo_km,'Topography',5,5, isLeft=True, latlon = False,  
#           xmin = 70., xmax = 105., ymin = 15., ymax = 45., prjstr = p_lcc)
