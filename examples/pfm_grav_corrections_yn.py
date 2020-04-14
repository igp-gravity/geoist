# -*- coding: utf-8 -*-
"""
Created on Sun Apr 12 22:46:16 2020

@author: chens
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
from geoist.pfm import normgra
from geoist import DATA_PATH 
# 1.读取数据
datapath = Path(Path(normgra.__file__).parent, 'data')
filename = Path(datapath, 'ynyx_grav.csv') #'ynyx_grav.csv') #ynp1_grav.csv
gradata = pd.read_csv(filename)
print('1. 重力数据已经读入，格式为：{}'.format(gradata.keys()))


# 2. 计算FGA
gradata['freeair'] = normgra.gamma_closed_form(gradata['lat'], gradata['elev'])
gradata['buglayer']  = normgra.bouguer_plate(gradata['elev'])
gradata['FGA']  = gradata['grav'] - gradata['freeair']
gradata['BGA_s']  = gradata['grav'] - gradata['freeair'] - gradata['buglayer'] 

# 2.1 正常场计算方法不同
#gradata['BGA_s1']  = gradata['grav'] - normgra.gamma_somigliana_free_air(gradata['lat'], gradata['elev']) - gradata['buglayer'] 

# 2.2 输出结果


print('2. 重力异常计算完成，已保存到：{}'.format(Path(DATA_PATH, 'ynp1_grav_anomaly.csv')))

# 3. 投影变换
import pyproj
p_jw = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
p_lcc = "+proj=lcc +lon_0=102.5 +lat_0=24.38 +lat_1=45 +ellps=WGS84 +datum=WGS84 +no_defs"

proj_xy = pyproj.Proj(p_lcc) #projection = pyproj.Proj(proj="merc", lat_ts=gradata['lat'].mean())
proj_coords = proj_xy(gradata['lon'].values, gradata['lat'].values)
gradata['x'] = proj_coords[0]
gradata['y'] = proj_coords[1]

proj_jw =  pyproj.Proj(p_jw) # 目标坐标系统
origin_lon, origin_lat = pyproj.transform(proj_xy, proj_jw, gradata['x'].values, gradata['y'].values)


fig, (ax0, ax1) = plt.subplots(nrows=1, ncols=2,figsize=(12,6))
ax0.set_title("Locations of gravity anomlay")
ax0.plot(gradata['lon'], gradata['lat'], "ok", markersize=1.5)
ax0.set_xlabel("Longitude")
ax0.set_ylabel("Latitude")
ax0.set_aspect("equal")

ax1.set_title("Projected coordinates of gravity anomlay")
ax1.plot(gradata['x'], gradata['y'], "ok", markersize=1.5)
ax1.set_xlabel("Easting (m)")
ax1.set_ylabel("Northing (m)")
ax1.set_aspect("equal")
plt.tight_layout()
plt.show()
# 4. 网格化
from geoist.gridder import spline, mask, trend

BGA = spline.Spline().fit(proj_coords, gradata['BGA_s'].values)

res_BGA = gradata['BGA_s'].values - BGA.predict(proj_coords) 

grid = BGA.grid(spacing=2e2, data_names=["BGAs"])
print('4. 网格化信息如下：',grid)
type(grid)

grid1 = mask.distance_mask(proj_coords, maxdist=5e2, grid=grid)
grid1.BGAs.plot.pcolormesh()

# 5. 去趋势

trend = trend.Trend(degree = 1).fit(proj_coords, gradata['BGA_s'].values)
print('4. 拟合的趋势系数：'.format(trend.coef_))

trend_values = trend.predict(proj_coords)
residuals = gradata['BGA_s'].values - trend_values

gradata['resBGA'] = residuals
gradata.to_csv(Path(DATA_PATH, 'ynp1_grav_anomaly.csv'), index = False, sep = ',')

ori = proj_coords #(origin_lon, origin_lat)
rBGA = spline.Spline().fit(ori, residuals)
grid2 = rBGA.grid(spacing=200, data_names=["resBGA"])


# grid3 = mask.distance_mask(ori, maxdist=5e2, grid=grid2)
plt.figure()
grid2.resBGA.plot.pcolormesh()
# grid3.resBGAs.plot.pcolormesh()

# 6. 数据转换
from geoist.others.utils import grid2Grid,map2DGrid
from geoist.pfm.grdio import grddata
rBGAg2d = grid2Grid(grid2.resBGA.easting, grid2.resBGA.northing, grid2.resBGA.values)
print(rBGAg2d)
g1out = grddata()
g1out.cols = rBGAg2d.getGeoDict().nx
g1out.rows = rBGAg2d.getGeoDict().ny
g1out.xmin = rBGAg2d.getGeoDict().xmin
g1out.xmax = rBGAg2d.getGeoDict().xmax
g1out.ymin = rBGAg2d.getGeoDict().ymin
g1out.ymax = rBGAg2d.getGeoDict().ymax
g1out.data0 = grid2.resBGA.values
g1out.export_surfer(Path(DATA_PATH, 'ynyx_bgar.grd'), False, 'ascii')

g1 = spline.Spline().fit(ori, gradata['FGA'])
gridx = g1.grid(spacing=200, data_names=["vals"])
g1out.data0 = gridx.vals.values
g1out.export_surfer(Path(DATA_PATH, 'ynyx_fga.grd'), False, 'ascii')

g2 = spline.Spline().fit(ori, gradata['BGA_s'])
gridx = g2.grid(spacing=200, data_names=["vals"])
g1out.data0 = gridx.vals.values
g1out.export_surfer(Path(DATA_PATH, 'ynyx_bgas.grd'), False, 'ascii')

g3= spline.Spline().fit(ori, gradata['elev'])
gridx = g3.grid(spacing=200, data_names=["vals"])
g1out.data0 = gridx.vals.values
g1out.export_surfer(Path(DATA_PATH, 'ynyx_elev.grd'), False, 'ascii')


#plt.figure()
#map2DGrid(None,rBGAg2d,'residual BGA', 200,200, isLeft=True)

plt.figure()
plt.imshow(grid2.resBGA.values[::-1])
np.save('d:\grid',grid2.resBGA.values[::-1])




