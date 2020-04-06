# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 20:13:29 2020

@author: chens
"""
import geoist as gi
import geoist.others.fetch_data as data
from geoist.others.fetch_data import _retrieve_file as downloadurl
from geoist.others.fetch_data import usgs_catalog
from geoist.catalog import Catalogue as cat
from geoist.catalog import QCmulti as cp
from geoist.catalog import QCreport as qc
from geoist.catalog import Smoothing as sm
from geoist.catalog import Selection as sel
from geoist.catalog import CatUtils as ct
from geoist.catalog import MapTools as mapt
from geoist.catalog import Exploration as exp
#import pdb

usgsfile = 'usgsca.csv'
localpath2 = usgs_catalog(usgsfile, '1970-01-01', '2020-01-01', '15','55','70','135',minmag = '5')
print(localpath2) 
dbusgs = cat.Database(usgsfile)
dbusgs.Import0(localpath2)
dbusgs.Info()

#平滑地震目录
#pdb.set_trace()
p = [(90.,20.),(90.,40.),(105.,40.),(105.,20.),(90.,20.)]
db3 = sel.AreaSelect(dbusgs,p)
P = ct.Polygon()
P.Load(p)

x1,y1,z1 = exp.GetHypocenter(db3)
wkt = ct.XYToWkt(P.x, P.y)
xsm, ysm, asm = sm.SmoothMFD(db3, 1., wkt, Delta=0.5)
cfg1 = {'Bounds': [90., 20., 105., 40.],
        'FigSize': [10., 12.],
        'Background': ['none',[0.9,0.8,0.6],[0.5,0.8,1.]],
        'Grid': [5., 5.]}
m1 = mapt.GeoMap(cfg1)
m1.BasePlot()
m1.MeshPlot(xsm, ysm, asm)
#m1.AreaPlot(P.x, P.y, Set=['y',0.5,'k',1])
#m1.PointPlot(xsm, ysm, Set=['o','b',2,1], Label='Grid')
m1.PointPlot(x1, y1, Set=['o','g',5,1], Label='全部')
m1.DrawGrid()
m1.Title('川滇地区地震目录高斯平滑')
m1.Show()

from geoist.catalog import Declusterer as declus
dbm, log1 = declus.WindowSearch(dbusgs, WinFun= declus.GardnerKnopoff, WinScale=1)
dbm.Info()


from geoist.catalog import MagRules as mr
sel.MagConvert(dbm,'*',['mb','Mb'],'Mw',mr.mb_Mw_Lin_DiGiacomo2015)
dbm.Info()

qc.pathname = gi.TEMP_PATH
# qc.network = 'usgs'
# qc.start_year = '2000'
# qc.end_year = '2020'
# qc.time_window = 2.0
# qc.dist_window = 15.0

# ## 地震目录质量检测
# pathname,prefix = qc.create_figures_new(qc.qcinit(), dbusgs)
# ## 生成HTML报告
# qc.generate_html(pathname,prefix,qc.to_show)

url = data.ispec_catalog_url
filename = '2020-03-25CENC-M4.7.dat'
localpath = downloadurl(url+filename, filename)
print(localpath) 
# catname = 'CENCM4.7'
# db2 = cat.Database(catname)
# header = ['Year', 'Month','Day','Hour','Minute','Second','Latitude','Longitude', 'Depth','MagType','MagSize','Log']
# db2.Import0(localpath, Header = header, Delimiter= ' ', flag = False)
# db2.SetField('LocCode','CENC')
# db2.SetField('MagCode','CENC4.7')
# #地震筛选
# lon = [70, 135]
# lat = [15, 55]
# db2.Filter('Latitude',lat[0],Opr='>=')
# db2.Filter('Latitude',lat[1],Opr='<=')
# db2.Filter('Longitude',lon[0],Opr='>=')
# db2.Filter('Longitude',lon[1],Opr='<=')
# db2.Info()

# outputname = cp.create_figures_new(db = [db2, dbusgs], pathname = gi.TEMP_PATH, 
#                                     startyear = 1970 , endyear = 2020, dhrs = 8)
# cp.generate_html(outputname, True)

