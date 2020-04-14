# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 08:28:28 2020

@author: chens
"""
from pathlib import Path
import pandas as pd
from geoist import EXAMPLES_PATH
from geoist import catalog
from geoist.gravity.gratools import Mapviz
from geoist.others.fetch_data import usgs_catalog
from geoist.catalog import Catalogue as cat


map1 = Mapviz(loc = 'cn')
map1.gmt_plot_base(map1.region)
 
usgsfile = 'usgscatx.csv' #下载保存到本地的文件名
localpath2 = usgs_catalog(usgsfile, '2009-01-01', '2014-01-02', map1.region[2],map1.region[3],map1.region[0],map1.region[1],minmag = '5')
print(localpath2)  #完整路径信息
data = pd.read_csv(localpath2)
map1.gmt_plot_xy(data)
filename = Path(EXAMPLES_PATH, 'mapviz1.png')
#map1.savefigure(filename)


#map1.fig.show()
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.io.shapereader import Reader
import matplotlib.pyplot as plt
import shapely.geometry as sgeom
from matplotlib.offsetbox import AnchoredText

datapath = Path(Path(catalog.__file__).parent,'data')



ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([72, 137, 10, 55])
ax.stock_img()
ax.coastlines()
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.OCEAN)
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.LAKES)
#ax.add_feature(cfeature.BORDERS, linestyle=':')
ax.add_feature(cfeature.RIVERS)

extent = (95, 108, 19.5, 44.5)
extent_box = sgeom.box(extent[0], extent[2], extent[1], extent[3])
xticks = range(72, 137, 5)
yticks = range(15, 55, 5)
ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True)

ax.add_geometries([extent_box], ccrs.PlateCarree(), color='white',
                  alpha=0.5, edgecolor='black', linewidth=2)

fname = Path(datapath, 'bou1_4l.shp')
f2name = Path(datapath, 'bou2_4l.shp')
faults = Path(datapath, 'gem_active_faults.shp')

ax.add_geometries(Reader(str(faults)).geometries(),
                  ccrs.PlateCarree(),facecolor = 'none',
                  edgecolor='red')

ax.add_geometries(Reader(str(f2name)).geometries(),
                  ccrs.PlateCarree(),  facecolor = 'none', 
                  edgecolor='gray', linestyle=':')

ax.add_geometries(Reader(str(fname)).geometries(),
                  ccrs.PlateCarree(),  facecolor = 'none', 
                  edgecolor='black')

scatter = ax.scatter(data.longitude, data.latitude,
           s= (0.2* 2 ** data.mag)**2, 
           c=data.depth / data.depth.max(), alpha=0.8,
           transform=ccrs.PlateCarree())


# produce a legend with the unique colors from the scatter
kw = dict(prop="colors", num= 6, fmt="{x:.0f} km",
          func=lambda s: s*data.depth.max())

legend1 = ax.legend(*scatter.legend_elements(**kw),
                    loc="lower left", title="Depth")
ax.add_artist(legend1)

kw = dict(prop="sizes", num=5, color=scatter.cmap(0.7), fmt="M {x:.1f}",
          func=lambda s: np.log2(np.sqrt(s)/0.2))
legend2 = ax.legend(*scatter.legend_elements(**kw),
                    loc="upper left", title="Mag")

SOURCE = 'GEOIST 2020'
LICENSE = 'MIT'
text = AnchoredText(r'$\mathcircled{{c}}$ {}; license: {}'
                    ''.format(SOURCE, LICENSE),
                    loc=4, prop={'size': 12}, frameon=True)
ax.add_artist(text)

plt.show()

# ax = plt.axes(projection=ccrs.Mollweide())
# ax.stock_img()
# plt.show()

# ax = plt.axes(projection=ccrs.PlateCarree())
# ax.stock_img()

# ny_lon, ny_lat = -75, 43
# delhi_lon, delhi_lat = 77.23, 28.61

# plt.plot([ny_lon, delhi_lon], [ny_lat, delhi_lat],
#          color='blue', linewidth=2, marker='o',
#          transform=ccrs.Geodetic(),
#          )

# plt.plot([ny_lon, delhi_lon], [ny_lat, delhi_lat],
#          color='gray', linestyle='--',
#          transform=ccrs.PlateCarree(),
#          )

# plt.text(ny_lon - 3, ny_lat - 12, 'New York',
#          horizontalalignment='right',
#          transform=ccrs.Geodetic())

# plt.text(delhi_lon + 3, delhi_lat - 12, 'Delhi',
#          horizontalalignment='left',
#          transform=ccrs.Geodetic())

# plt.show()