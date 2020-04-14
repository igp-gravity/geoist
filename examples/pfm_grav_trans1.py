# -*- coding: utf-8 -*-
"""
Created on Sun Apr 12 22:46:16 2020

@author: chens
"""

import matplotlib.pyplot as plt
from pathlib import Path
from geoist import DATA_PATH 
from geoist.others.gdal import GDALGrid
from geoist.others.utils import map2DGrid, Grid2Xyz, xyz2Grid, grid2srf

# 1.读取数据

filename1 = Path(DATA_PATH, 'ynyx_elev.grd')
gd1, ff = GDALGrid.getFileGeoDict(filename1)
grid1 = GDALGrid.load(filename1,gd1)
filename2 = Path(DATA_PATH, 'ynyx_fga.grd')
gd2, ff = GDALGrid.getFileGeoDict(filename2)
grid2 = GDALGrid.load(filename2,gd2)
filename3 = Path(DATA_PATH, 'ynyx_bgas.grd')
gd3, ff = GDALGrid.getFileGeoDict(filename3)
grid3 = GDALGrid.load(filename3,gd3)
filename4 = Path(DATA_PATH, 'ynyx_bgar.grd')
gd4, ff = GDALGrid.getFileGeoDict(filename4)
grid4 = GDALGrid.load(filename4,gd4)



from geoist.pfm import pftrans
from geoist.vis import giplt

x, y ,elev = Grid2Xyz(grid1)
x, y ,fga = Grid2Xyz(grid2)
x, y ,bgas = Grid2Xyz(grid3)
x, y ,bgar = Grid2Xyz(grid4)
shape = (grid3.getGeoDict().nx, grid3.getGeoDict().ny)
height = 1000  # How much higher to go
bgas_contf = pftrans.upcontinue(x, y, bgas, shape, height)

bgas_dx = pftrans.derivx(x, y, bgas_contf, shape)
bgas_dy = pftrans.derivy(x, y, bgas_contf, shape)
bgas_dz = pftrans.derivz(x, y, bgas_contf, shape)
bgas_tga = pftrans.tga(x, y, bgas_contf, shape)

import numpy as np
from geoist.others.geodict import GeoDict

xmin = np.min(x)
xmax = np.max(x)
ymin = np.min(y)
ymax = np.max(y)
nx = len(set(x))
ny = len(set(y))
dx = (xmax-xmin)/(nx-1)
dy = (ymax-ymin)/(ny-1)      
geodict = {'xmin':xmin,'xmax':xmax,'ymin':ymin,'ymax':ymax,'dx':dx,'dy':dy,'nx':nx,'ny':ny}
    
bgascgrd = xyz2Grid(x, y, bgas_contf)
grid2srf(bgascgrd)

args = dict(shape=shape, levels=20, cmap=plt.cm.RdBu_r)
fig, axes = plt.subplots(2, 2, figsize=(10, 10))
plt.sca(axes[0, 0])
plt.title("Derivative of BGA in X")
plt.axis('scaled')
giplt.contourf(x, y, bgas_dx, **args)
plt.colorbar(pad=0).set_label('mGal/m')
giplt.m2km()

plt.sca(axes[0, 1])
plt.title("Derivative of BGA in Y")
plt.axis('scaled')
giplt.contourf(x, y, bgas_dy, **args)
plt.colorbar(pad=0).set_label('mGal/m')
giplt.m2km()

plt.sca(axes[1, 0])
plt.title("Derivative of BGA in Z")
plt.axis('scaled')
giplt.contourf(x, y, bgas_dz, **args)
plt.colorbar(pad=0).set_label('mGal/m')
giplt.m2km()

plt.sca(axes[1, 1])
plt.title("Total gradient amplitude of BGA")
plt.axis('scaled')
giplt.contourf(x, y, bgas_tga, **args)
plt.colorbar(pad=0).set_label('mGal/m')
giplt.m2km()