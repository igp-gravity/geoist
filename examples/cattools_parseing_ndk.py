"""
EXAMPLE 9 - Parseing GCMT bulletin (Ndk format)
"""
from os.path import dirname
from geoist.cattools import Parsers as Par
from geoist.cattools import Exploration as Exp
from geoist.cattools import MapTools as Map
from geoist.cattools import Selection as Sel
#-----------------------------------------------------------------------------------------
# Import catalogue
pathname = dirname(__file__)
Db = Par.Database('GCMT')
Db.ImportNdk(pathname+'/data/jan76_dec13.ndk')

#-----------------------------------------------------------------------------------------
# Search Area (Africa)
p = [(-20.,30.),(60.,30.),(60.,60.),(-20.,60.)]

Db2 = Sel.AreaSelect(Db,p)

#-----------------------------------------------------------------------------------------
# Plot map

cfg = {'Bounds': [-20., 30., 60., 60.],
       'FigSize': [8., 6.],
       'Background': ['none',[0.9,0.8,0.6],[0.5,0.8,1.]],
       'Grid': [10., 10.]}

M = Map.GeoMap(cfg)

M.BasePlot()
M.DrawBounds()
M.DrawGrid()

MagTab = [[3,4,['d','w',2,1]],
          [4,5,['^','y',4,1]],
          [5,6,['o','g',6,1]],
          [6,7,['s','r',8,1]],
          [7,8,['p','m',10,1]]]

for mt in MagTab:

  # Selection by magnitude range
  DbC = Sel.MagRangeSelect(Db, mt[0], mt[1])
  x,y,z = Exp.GetHypocenter(DbC)
  M.PointPlot(x, y, Set=mt[2])

M.Title('Example - cattools parseing ndk')
M.Show()

#M.SaveFig(pathname+'/data/example9.png')
