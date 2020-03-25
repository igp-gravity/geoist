"""
EXAMPLE 7 - Prime analysis
"""
from os.path import dirname
from geoist.catalog import Parsers as Par
from geoist.catalog import Exploration as Exp
from geoist.catalog import MapTools as Map
from geoist.catalog import Selection as Sel
#-----------------------------------------------------------------------------------------
# Import catalogue
pathname = dirname(__file__)
Db = Par.Database('ISC-HOM')
Db.ImportIsf(pathname+'/data/isc-rev-africa.isf')

#-----------------------------------------------------------------------------------------
# Split primes

DbP, DbN = Sel.SplitPrime(Db)

xn = DbN.Extract('Longitude')
yn = DbN.Extract('Latitude')

xp = DbP.Extract('Longitude')
yp = DbP.Extract('Latitude')

#-----------------------------------------------------------------------------------------
# Plot map

cfg = {'Bounds': [10., -40., 60., 20.],
       'FigSize': [8., 6.],
       'Background': ['none',[0.9,0.8,0.6],[0.5,0.8,1.]],
       'Grid': [10., 10.]}

M = Map.GeoMap(cfg)

M.BasePlot()
M.DrawBounds()
M.DrawGrid()

M.PointPlot(xp, yp, Set=['o','b',5,1], Label='Primes')
M.PointPlot(xn, yn, Set=['o','r',5,1], Label='Not Primes')

M.Legend()
M.Title('Example - prime analysis')
M.Show()

#M.SaveFig(pathname+'/data/example7.png')
