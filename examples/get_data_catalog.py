# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 17:11:22 2020

@author: chens
"""

#from os.path import dirname
import numpy as np
import geoist as gi
import geoist.others.fetch_data as data
from geoist.others.fetch_data import _retrieve_file as downloadurl

from geoist.catalog import QCreport as qc
from geoist.catalog import QCmulti as cp
from geoist.catalog import Catalogue as cat
from geoist.catalog import Exploration as exp
from geoist.catalog import MapTools as mapt
from geoist.catalog import Selection as sel
from geoist.catalog import Seismicity as sem
from geoist.catalog import Declusterer as declus
from geoist.catalog import Smoothing as sm
from geoist.catalog import CatUtils as ct

## 下载地震目录
# pathname = dirname(__file__)
# print(pathname)
url = data.ispec_catalog_url
print(url)
filename = '2020-03-25CENC-M4.7.dat'

localpath = downloadurl(url+filename, filename)
print(localpath) #文件路径

## 建立地震目录数据库
catname = 'CENCM4.7'
db2 = cat.Database(catname)
header = ['Year', 'Month','Day','Hour','Minute','Second','Latitude','Longitude', 'Depth','MagType','MagSize','Log']
db2.Import0(localpath, Header = header, Delimiter= ' ', flag = False)
db2.Info()


db2.SetField('LocCode','CENC')
db2.SetField('MagCode','CENC4.7')


#地震筛选
# Search Area (China) using internal filter
lon = [70, 135]
lat = [15, 55]
db2.Filter('Latitude',lat[0],Opr='>=')
db2.Filter('Latitude',lat[1],Opr='<=')
db2.Filter('Longitude',lon[0],Opr='>=')
db2.Filter('Longitude',lon[1],Opr='<=')

db2.Info()
#二维时间序列图
exp.MagTimePlot(db2)
exp.MagTimeBars(db2)
exp.RateDensityPlot(db2)
# G-R关系
enum, mbin =exp.GetKeyHisto(db2,'MagSize',Bnum=20, Norm=False) 
minc= (max(mbin)-min(mbin))/10.
#拟合b值
a,b = sem.MfdOptimize(enum, mbin, minc, max(mbin))  
print('b-value=',b) 
#复发概率
sem.MfdPlot(a,b, max(mbin),Enum=enum, Ecum=np.cumsum(enum[::-1])[::-1], Mbin=mbin, Minc=[minc])

## 去余震
dbm, log1 = declus.WindowSearch(db2, WinFun= declus.GardnerKnopoff, WinScale=1)
dbm.Info()
## 震中分布图
x1,y1,z1 = exp.GetHypocenter(db2)
x2,y2,z2 = exp.GetHypocenter(dbm)

cfg = {'Bounds': [70., 15., 135., 55.],
        'FigSize': [16., 12.],
        'Background': ['none',[0.9,0.8,0.6],[0.5,0.8,1.]],
        'Grid': [10., 10.]}

M = mapt.GeoMap(cfg)
M.BasePlot()
M.DrawBounds()
M.DrawGrid()
#震中分布图
M.PointPlot(x1, y1, Set=['o','g',5,1], Label='全部')
M.PointPlot(x2, y2, Set=['*','r',2,1], Label='去余震')
M.Legend()
M.Title('中国及邻区震中分布图')
M.Show()

#平滑地震目录
p = [(90.,20.),(90.,40.),(105.,40.),(105.,20.),(90.,20.)]
db3 = sel.AreaSelect(db2,p)
P = ct.Polygon()
P.Load(p)

db3.Info()
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

## 得到系统路径和记录日志
#print(gi.EXAMPLES_PATH, gi.DATA_PATH, gi.TEMP_PATH)
nwcat = qc.pathname+'\\mwcat1900utc.csv'
print(qc.pathname+'\\cn-cat-mw.txt')
qc.pathname = gi.TEMP_PATH
qc.network = catname
qc.start_year = '1970'
qc.end_year = '2020'
qc.time_window = 2.0
qc.dist_window = 15.0
##gi.__verison__
gi.log.info(catname+'/catalog qctest')

## 地震目录质量检测
pathname,prefix = qc.create_figures_new(qc.qcinit(), db2)

## 生成HTML报告
qc.generate_html(pathname,prefix,qc.to_show)

## 地震目录对比
dbm.Header['Name'] = 'CENCm'
db2.Header['Name'] = 'CENC4.7'

## 建立地震目录数据库
catname = 'cnmw'
localpath = nwcat
db6 = cat.Database(catname)
header = ['Year', 'Month','Day','Hour','Minute','Second','Latitude','Longitude','MagSize','Depth','Log']
db6.Import0(localpath, Header = header, Delimiter= ',', flag = False)

db6.SetField('MagType', 'Mw')
db6.Info()
outputname = cp.create_figures_new(db = [db2, db6], pathname = gi.TEMP_PATH, 
                                   startyear = 1970 , endyear = 2015, dhrs = 8)

## 生成HTML报告
no_output_matches = True
cp.generate_html(outputname, no_output_matches)
