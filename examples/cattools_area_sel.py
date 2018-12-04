"""
EXAMPLE 6 - Area selection
"""

from os.path import dirname

from geoist.cattools import Catalogue as Cat
from geoist.cattools import Exploration as Exp
from geoist.cattools import Selection as Sel
from geoist.cattools import MapTools as Map
from geoist.cattools import CatUtils as Ct

pathname = dirname(__file__)
print(pathname+"/data/isc-rev-africa-select.bin")
#-----------------------------------------------------------------------------------------
# Import catalogue
Db1 = Cat.Database()
Db1.Load(pathname+'/data/isc-rev-africa-select.bin')

#-----------------------------------------------------------------------------------------
# Filtering by polygon

# As WKT Format or List
CASE = 1


if CASE == 1:
  p = [(20.,-20.),(40.,-20.),(40.,0.),(20.,0.)]
  Db2 = Sel.AreaSelect(Db1,p)
  P = Ct.Polygon()
  P.Load(p)

if CASE == 2:
  p = pathname+"/data/area.xy"
  Db2 = Sel.AreaSelect(Db1,p,File='xy')
  P = Ct.Polygon()
  P.Import(p)

if CASE == 3:
  p = pathname+"/data/area.wkt"
  Db2 = Sel.AreaSelect(Db1,p,File='wkt')
  P = Ct.Polygon()
  P.Import(p,Type='wkt')

#-----------------------------------------------------------------------------------------
# Get Coordinates

x1,y1,z1 = Exp.GetHypocenter(Db1)
x2,y2,z2 = Exp.GetHypocenter(Db2)

#-----------------------------------------------------------------------------------------
# Map Plot

# Reshape polygon (to plot area)
P.x.append(P.x[0])
P.y.append(P.y[0])

# Plot map
cfg = {'Bounds': [10., -40., 60., 20.],
       'FigSize': [8., 6.],
       'Background': ['none',[0.9,0.8,0.6],[0.5,0.8,1.]],
       'Grid': [10., 10.]}

M = Map.GeoMap(cfg)

M.BasePlot()
M.DrawBounds()
M.DrawGrid()

M.PointPlot(x1, y1, Set=['o','g',5,1], Label='All')
M.PointPlot(x2, y2, Set=['o','r',5,1], Label='Selection')

M.AreaPlot(P.x, P.y, Set=['y',0.5,'k',1])

M.Legend()
M.Title('Example-cattools_area_sel')
M.Show()

#M.SaveFig(pathname+"/data/cattools_area_sel.jpg")
