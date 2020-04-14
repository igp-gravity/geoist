# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 13:40:07 2020

@author: chens
"""

import matplotlib.pyplot as plt
import numpy as np
from geoist.others.grid2d import Grid2D
from geoist.others.geodict import GeoDict
from mpl_toolkits.basemap import Basemap
import warnings
warnings.simplefilter("ignore")

def map2DGrid(ax, grid, tstr, xlen=1.0, ylen=1.0, isLeft=False):
    """
    grid is a Grid2D object 
    """
    xmin,xmax,ymin,ymax = grid.getBounds()
    pdata = grid.getData()
    nr,nc = pdata.shape
    lonrange = np.linspace(xmin,xmax,num=nc)
    latrange = np.linspace(ymin,ymax,num=nr)
    lon,lat = np.meshgrid(lonrange,latrange)
    latmean = np.mean([ymin,ymax])
    lonmean = np.mean([xmin,xmax])
    m = Basemap(llcrnrlon=xmin,llcrnrlat=ymin,urcrnrlon=xmax,urcrnrlat=ymax,\
            rsphere=(6378137.00,6356752.3142),\
            resolution='c',area_thresh=1000.,projection='lcc',\
            lat_1=latmean,lon_0=lonmean,ax=ax)
    # draw coastlines and political boundaries.
    m.drawcoastlines()
    #m.drawcountries()
    #m.drawstates()
    lons = np.arange(xmin,xmax,xlen)
    lats = np.arange(ymin,ymax,ylen)
    if isLeft:
        labels = labels=[1,0,0,0]
    else:
        labels = labels=[0,0,0,0]
    m.drawparallels(lats,labels=labels,color='white',fmt='%.1f') # draw parallels
    m.drawmeridians(lons,labels=[0,0,0,1],color='white',fmt='%.1f') # draw meridians
    pmesh = m.pcolormesh(lon,lat,np.flipud(grid.getData()),latlon=True)
    #plt.hold(True)
    ax.set_title(tstr)
    m.colorbar(pmesh)
    
    
xmin = 118.5
xmax = 120.5
ymin = 32.0
ymax = 34.0
xdim = 0.25
ydim = 0.25

ncols = len(np.arange(xmin,xmax+xdim,xdim))
nrows = len(np.arange(ymin,ymax+ydim,ydim))
data = np.arange(0,nrows*ncols)
data.shape = (nrows,ncols)

geodict = {'xmin':xmin,
           'xmax':xmax,
           'ymin':ymin,
           'ymax':ymax,
           'dx':xdim,
           'dy':ydim,
           'nx':ncols,
           'ny':nrows,}

grid = Grid2D(data,GeoDict(geodict))

plt.figure()
plt.imshow(grid.getData(),interpolation='nearest')
plt.colorbar()

lat,lon = grid.getLatLon(4,4)
print('The coordinates at the center of the grid are %.3f,%.3f' % (lat,lon))
value = grid.getValue(lat,lon)
print('The data value at the center of the grid is %i' % (value))

xmin = 119.2
xmax = 119.8
ymin = 32.7
ymax = 33.3
dx = 0.33
dy = 0.33
nx = len(np.arange(xmin,xmax+xdim,xdim))
ny = len(np.arange(ymin,ymax+ydim,ydim))
sdict = {'ymin':32.7,'ymax':33.3,'xmin':119.2,'xmax':119.8,'dx':0.33,'dy':0.33,'nx':nx,'ny':ny}
grid2 = grid.interpolateToGrid(GeoDict(sdict))

plt.figure()
plt.imshow(grid2.getData(),interpolation='nearest'); #'linear', 'cubic', 'nearest'
plt.colorbar()

# from geoist.others.gmt import GMTGrid
# gd0, ff = GMTGrid.getFileGeoDict('d:\\out.nc')
# grid0 = GMTGrid.load('d:\\out.nc',gd0)
# plt.figure()
# plt.imshow(grid0.getData())
# plt.colorbar()

from geoist.others.gdal import GDALGrid
# gd1, ff = GDALGrid.getFileGeoDict('d:\\out-arc.grd')
# grid = GDALGrid.load('d:\\out-arc.grd',gd1)
# plt.figure()
# plt.imshow(grid.getData())
# plt.colorbar()

filename1 = 'D:\\research-work\\Topex\\WGM2012\\WGM2012_ETOPO1_ponc_2min.grd'
filename2 = 'D:\\research-work\\Topex\\WGM2012\\WGM2012_Freeair_ponc_2min.grd'
filename3 = 'D:\\research-work\\Topex\\WGM2012\\WGM2012_Bouguer_ponc_2min.grd'
filename4 = 'D:\\research-work\\Topex\\WGM2012\\WGM2012_Isostatic_ponc_2min.grd'


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
topo = GDALGrid.load(filename1, samplegeodict = gd1, resample = True)
fag = GDALGrid.load(filename2, samplegeodict = gd1, resample = True)
bug = GDALGrid.load(filename3, samplegeodict = gd1, resample = True)
iso = GDALGrid.load(filename4, samplegeodict = gd1, resample = True)
plt.figure()
plt.imshow(topo.getData())
plt.colorbar()

#sdict = {'ymin':15.0,'ymax':45.0,'xmin':70.0,'xmax':105.0,'dx':0.5,'dy':0.5,'nx':800,'ny':600}
#topo2 = topo.interpolateToGrid(GeoDict(sdict, adjust='res'))

# fag = GDALGrid.load(filename2,gd2)
# bug = GDALGrid.load(filename3,gd3)



# # 投影，取数据，内插，外插，切片

fig, axes = plt.subplots(nrows=2,ncols=2,figsize=(12,12))
fig.tight_layout()
map2DGrid(axes[0,0],topo,'Topography',10,10, isLeft=True)
map2DGrid(axes[0,1],fag,'Free-air gravity',10,10)
map2DGrid(axes[1,0],bug,'Bouguer gravity',10,10, isLeft=True)
map2DGrid(axes[1,1],iso,'Isostatic gravity',10,10)
