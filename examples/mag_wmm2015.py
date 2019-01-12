# -*- coding: utf-8 -*-
"""
Created on Thu Jan 10 18:34:07 2019
计算WMM2015模型，WMM.cof文件需要放到与py相同目录
@author: chens
"""

import numpy as np
from pathlib import Path
import xarray
import ctypes as ct
import sys
import datetime
from matplotlib.pyplot import figure

#libwmm = ct.cdll.LoadLibrary(str('D:\\MyWorks\\WMM2015-master\\wmm15.dll'))
libwmm = ct.cdll.LoadLibrary(str('D:\\MyWorks\\WMM2015-master\\noaa.dll'))

def noaa(glats: np.ndarray, glons: np.ndarray, alt_km: float, yeardec: float, mod = 'wmm') -> xarray.Dataset:

    glats = np.atleast_2d(glats).astype(float)  # to coerce all else to float64
    glons = np.atleast_2d(glons)

    assert glats.shape == glons.shape

    mag = xarray.Dataset(coords={'glat': glats[:, 0], 'glon': glons[0, :]})
    north = np.empty(glats.size)
    east = np.empty(glats.size)
    down = np.empty(glats.size)
    total = np.empty(glats.size)
    decl = np.empty(glats.size)
    incl = np.empty(glats.size)

    for i, (glat, glon) in enumerate(zip(glats.ravel(), glons.ravel())):

        x = ct.c_double()
        y = ct.c_double()
        z = ct.c_double()
        T = ct.c_double()
        D = ct.c_double()
        mI = ct.c_double()
        
        if mod == 'wmm':
            ret = libwmm.wmmsub(ct.c_double(glat),
                            ct.c_double(glon),
                            ct.c_double(alt_km),
                            ct.c_double(yeardec),
                            ct.byref(x), ct.byref(y), ct.byref(z),
                            ct.byref(T), ct.byref(D), ct.byref(mI))      
        else:
            ret = libwmm.emmsub(ct.c_double(glat),
                            ct.c_double(glon),
                            ct.c_double(alt_km),
                            ct.c_double(yeardec),
                            ct.byref(x), ct.byref(y), ct.byref(z),
                            ct.byref(T), ct.byref(D), ct.byref(mI))
        #print(ret)
        assert ret == 0

        north[i] = x.value
        east[i] = y.value
        down[i] = z.value
        total[i] = T.value
        decl[i] = D.value
        incl[i] = mI.value

    mag['north'] = (('glat', 'glon'), north.reshape(glats.shape))
    mag['east'] = (('glat', 'glon'), east.reshape(glats.shape))
    mag['down'] = (('glat', 'glon'), down.reshape(glats.shape))
    mag['total'] = (('glat', 'glon'), total.reshape(glats.shape))
    mag['incl'] = (('glat', 'glon'), incl.reshape(glats.shape))
    mag['decl'] = (('glat', 'glon'), decl.reshape(glats.shape))

    mag.attrs['time'] = yeardec

    return mag

def plotwmm(mag: xarray.Dataset):

    fg = figure()
    ax = fg.subplots(1, 2, sharey=True)
    fg.suptitle('WMM2015  {}'.format(mag.time))
    h = ax[0].contour(mag.glon, mag.glat, mag.decl, range(-90, 90+20, 20))
    ax[0].clabel(h, inline=True, fmt='%0.1f')
    ax[0].set_title('Magnetic Declination [degrees]')

    h = ax[1].contour(mag.glon, mag.glat, mag.incl, range(-90, 90+20, 20))
    ax[1].clabel(h, inline=True, fmt='%0.1f')
    ax[1].set_title('Magnetic Inclination [degrees]')

    ax[0].set_ylabel('Geographic latitude (deg)')
    for a in ax:
        a.set_xlabel('Geographic longitude (deg)')


from geoist.others.scidates import datetime2yeardec
dt = datetime.datetime(2012, 7, 12, 12)
print(datetime2yeardec(dt))
mag = noaa(45.5, 105.6, 0.2, datetime2yeardec(dt), mod='emm')
#print(mag.north.item())
#print(mag.east.item())
#print(mag.down.item())
print("F:",mag.total.item()) #F
print("D:",mag.decl.item())  #D
print("I:",mag.incl.item())  #I

from matplotlib.pyplot import show
lon, lat = np.meshgrid(np.arange(-180, 180+10, 10), np.arange(-90, 90+10, 10))
mag = noaa(lat, lon, 0, 2015)
plotwmm(mag)
show()