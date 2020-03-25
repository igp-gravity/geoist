"""
 Earthquake Catalog Analysis
"""
from os.path import dirname
import numpy as np
#local import
from geoist.catalog import Catalogue as Cat
from geoist.catalog import Exploration as Exp
from geoist.catalog import MapTools as Map
from geoist.catalog import Selection as Sel
from geoist.catalog import Seismicity as Sem
from geoist.catalog import Declusterer as Declus
from geoist.catalog import Smoothing as Sm
from geoist.catalog import CatUtils as Ct
#-----------------------------------------------------------------------------------------
pathname = dirname(__file__)
H = ['Id','','Year','Month','Day','Hour','Minute','Second',
     'Longitude','Latitude','','','','Depth','DepError',
     'MagSize','MagError','','','','','','','','','']

Db = Cat.Database('ISC-GEM')
Db.Import(pathname+'/data/isc-gem-v3.csv',Header=H, SkipLine=1, Delimiter=',')

Db.SetField('LocCode','ISC-GEM')
Db.SetField('MagCode','ISC-GEM')
Db.SetField('MagType','MW')
#-----------------------------------------------------------------------------------------
# Search Area (China) using internal filter
lon = [70, 135]
lat = [15, 55]
#地震筛选
Db.Filter('Latitude',lat[0],Opr='>=')
Db.Filter('Latitude',lat[1],Opr='<=')
Db.Filter('Longitude',lon[0],Opr='>=')
Db.Filter('Longitude',lon[1],Opr='<=')
Exp.AgencyReport(Db, 'L')
#二维时间序列图
Exp.MagTimePlot(Db)
Exp.MagTimeBars(Db)
Exp.RateDensityPlot(Db)
# G-R关系
Enum, Mbin =Exp.GetKeyHisto(Db,'MagSize',Bnum=10, Norm=False) 
Minc= (max(Mbin)-min(Mbin))/10.
#拟合b值
a,b = Sem.MfdOptimize(Enum, Mbin, Minc, max(Mbin))  
print('b-value=',b) 
#复发概率
Sem.MfdPlot(a,b, max(Mbin),Enum=Enum, Ecum=np.cumsum(Enum[::-1])[::-1], Mbin=Mbin, Minc=[Minc])
#重复事件监测
Log = Sel.MergeDuplicate(Db,Twin=60.,Swin=50.,Log=1)
Exp.DuplicateCheck(Log)
#去余震
Dbm, Log1 = Declus.WindowSearch(Db)
#目录摘要
Dbm.Info()
#震中提取
x1,y1,z1 = Exp.GetHypocenter(Db)
x2,y2,z2 = Exp.GetHypocenter(Dbm)

p = [(90.,20.),(90.,40.),(105.,40.),(105.,20.),(90.,20.)]
P = Ct.Polygon()
P.Load(p)

cfg = {'Bounds': [70., 15., 135., 55.],
       'FigSize': [8., 6.],
       'Background': ['none',[0.9,0.8,0.6],[0.5,0.8,1.]],
       'Grid': [10., 10.]}

M = Map.GeoMap(cfg)
M.BasePlot()
M.DrawBounds()
M.DrawGrid()
#震中分布图
M.PointPlot(x1, y1, Set=['o','g',5,1], Label='All')
M.PointPlot(x2, y2, Set=['*','r',2,1], Label='Main')
M.AreaPlot(P.x, P.y, Set=['y',0.5,'k',1])
#平滑地震目录
wkt = Ct.XYToWkt(P.x, P.y)
xsm, ysm, asm = Sm.SmoothMFD(Db, 1., wkt, Delta=0.5)
#M.PointPlot(xsm, ysm, Set=['o','b',2,1], Label='Grid')
M.MeshPlot(xsm, ysm, asm)
M.Legend()
M.Title('Earthquakes in China')
M.Show()

#print('dump to:'+pathname+'/data/isc-gem-v3.bin')
#Db.Dump(pathname+'/data/isc-gem-v3.bin')
