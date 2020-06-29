#!/usr/bin/env python
# -*- coding: utf-8 -*-


from shapely import wkt, geometry
import scipy as scp
import numpy as np
import math as ma
import re

#-----------------------------------------------------------------------------------------

def LocationInit():

  L = {'Year': None,
       'Month': None,
       'Day': None,
       'Hour': None,
       'Minute': None,
       'Second': None,
       'Latitude': None,
       'Longitude': None,
       'Depth': None,
       'SecError': None,
       'LatError': None,
       'LonError': None,
       'DepError': None,
       'LocCode': None,
       'Prime': False}

  return L

#-----------------------------------------------------------------------------------------

def MagnitudeInit():

  M = {'MagSize': None,
       'MagError': None,
       'MagType': None,
       'MagCode': None}

  return M

#-----------------------------------------------------------------------------------------

def CastValue(key, value):

  C = {'Id': 's',
       'Log': 's',
       'Year': 'i',
       'Month': 'i',
       'Day': 'i',
       'Hour': 'i',
       'Minute': 'i',
       'Second': 'f',
       'Latitude': 'f',
       'Longitude': 'f',
       'Depth': 'f',
       'SecError': 'f',
       'LatError': 'f',
       'LonError': 'f',
       'DepError': 'f',
       'LocCode': 's',
       'Prime': 'b',
       'MagSize': 'f',
       'MagError': 'f',
       'MagType': 's',
       'MagCode': 's'}

  if not IsEmpty(value):
    if C[key] == 'i':
      value = int(value)
    if C[key] == 'f':
      value = float(value)
    if C[key] == 's':
      value = str(value)
    if C[key] == 'b':
      value = bool(value)
  else:
    value = None

  return value

#-----------------------------------------------------------------------------------------

def KeyGroup(key):

  L = LocationInit()
  M = MagnitudeInit()
  G = []

  if key in L.keys():
    G = 'Location'

  if key in M.keys():
    G = 'Magnitude'

  return G

#-----------------------------------------------------------------------------------------

def IsEmpty(number):

  C0 = (number == [])
  C1 = (number == '')
  C2 = (number != number)
  C3 = (number == None)
  C4 = (number == 'None')

  return (C0 or C1 or C2 or C3 or C4)

#-----------------------------------------------------------------------------------------

def IsType(value, dtype):

  Out = False

  if dtype in ['Bool','bool','B','b']:
    if type(value) == bool:
      Out = True
  if dtype in ['Int','int','I','i']:
    if type(value) == int:
      Out = True
  if dtype in ['Float','float','F','f']:
    if type(value) == float:
      Out = True
  if dtype in ['String','string','S','s']:
    if type(value) == str:
      Out = True
  if dtype in ['List','list','L','l']:
    if type(value) == list:
      Out = True
  if dtype in ['Tuple','tuple','T','t']:
    if type(value) == tuple:
      Out = True

  return Out

#-----------------------------------------------------------------------------------------

def WgsDistance(Lat1, Lon1, Lat2, Lon2):
  """Author: Salvador Dali
  http://stackoverflow.com/users/1090562/salvador-dali
  
  """

  p = 0.017453292519943295

  c1 = ma.cos((Lat2 - Lat1) * p)
  c2 = ma.cos(Lat1 * p)
  c3 = ma.cos(Lat2 * p)
  c4 = ma.cos((Lon2 - Lon1) * p)

  a = 0.5 - c1/2 + c2 * c3 * (1 - c4) / 2

  return 12742 * ma.asin(ma.sqrt(a))

#-----------------------------------------------------------------------------------------

def WgsToXY (Lat, Lon, Km=True):
  """Approximate conversion using sinusoidal projection.

  """

  earth_radius = 6371009. # in meters
  y_dist = np.pi * earth_radius / 180.0

  Lat = np.array(Lat)
  Lon = np.array(Lon)

  y = Lat * y_dist
  x = Lon * y_dist * np.cos(np.radians(Lat))

  # Converting to Km
  if Km:
    y /= 1000.
    x /= 1000.

  return x, y

#-----------------------------------------------------------------------------------------

def ConcaveHull (X, Y):

  points = zip(X, Y)
  hull = scp.spatial.ConvexHull(points)
  ChX = points[hull.vertices,0]
  ChY = points[hull.vertices,1]

  return ChX, ChY

#-----------------------------------------------------------------------------------------

def SphericalMesh(Delta, Km=False):
  """
  Produce a shperical mesh using golder spiral algorithm.
  Distance is the average beween nearby points (degree by default).

  """

  if Km:
    # Distance is in Km
    Rad = 6371. # Approximated Earth radius
    N = np.rint((4*np.pi*Rad**2)/(Delta**2))
  else:
    # Distance is in Degree
    N = (4*np.pi)/(np.deg2rad(Delta)**2)

  I = np.arange(0, N, dtype=float) + 0.5

  Phi = np.arccos(1 - 2*I/N)
  Theta = np.pi * (1 + 5**0.5) * I

  # Conversion to Lat/Lon
  Lat = np.rad2deg(Phi) - 90.
  Lon = np.rad2deg(Unwrap(Theta))

  return Lon, Lat

#-----------------------------------------------------------------------------------------

def Unwrap(Angle):
  """
  Unwrap phase angle.
  Note: Angle must be a numpy array

  """

  return Angle-(2.*np.pi)*((Angle+np.pi)//(2*np.pi))

#-----------------------------------------------------------------------------------------

def LeapCheck(Year):

  C0 = (Year % 4 == 0)
  C1 = (Year % 100 != 0)
  C2 = (Year % 400 == 0)

  return (C0 and C1) or C2

def LeapNum(Year):

  N0 = (Year-1)//4
  N1 = (Year-1)//100
  N2 = (Year-1)//400

  return N0 - N1 + N2

#-----------------------------------------------------------------------------------------

def DateToSec(Year, Month, Day, Hour, Minute, Second):

  if Year < 1:
    print('Warning: Year must be > 1')
    return None

  if not Year: Year = 1.
  if not Month: Month = 1.
  if not Day: Day = 1.
  if not Hour: Hour = 0.
  if not Minute: Minute = 0.
  if not Second: Second = 0.

  DSEC = 24.*3600.
  YDAYS = 365.

  if LeapCheck(Year):
    MDAYS = [0.,31.,60.,91.,121.,152.,182.,213.,244.,274.,305.,335.]
  else:
    MDAYS = [0.,31.,59.,90.,120.,151.,181.,212.,243.,273.,304.,334.]

  YSec = (Year-1)*YDAYS*DSEC
  YSec += LeapNum(Year)*DSEC

  MSec = MDAYS[int(Month)-1]*DSEC

  DSec = (Day-1)*DSEC

  Sec = YSec + MSec + DSec + Hour*3600.+ Minute*60. + Second*1.

  return Sec

#-----------------------------------------------------------------------------------------

def WktToXY(WktString):
  """
  NOTE:
  1) Internal polygons are not considered
  2) Multi objects are simply concatenated

  """

  from array import array

  def multicoords(WktObj):
    """
    Subfunction to iterate though multi objects

    """
    X = array('d', [])
    Y = array('d', [])
    for M in WktObj:
      try:
        x,y = M.coords.xy
      except:
        x,y = M.exterior.coords.xy
      X += x
      Y += y
    return X, Y

  # Loading WKT
  WktObj = wkt.loads(WktString)

  if WktObj.type in ['Point', 'LineString', 'Polygon']:
    try:
      X,Y = WktObj.coords.xy
    except:
      X,Y = WktObj.exterior.coords.xy

  if WktObj.type in ['MultiPoint', 'MultiLineString', 'MultiPolygon']:
    X,Y = multicoords(WktObj)

  # Conversion to list
  X = [x for x in X]
  Y = [y for y in Y]

  return X, Y

#---------------------------------------------------------------------------------------

def XYToWkt(X, Y):

  WktString = 'POLYGON (('

  Data = []
  for x, y in zip(X, Y):
    Data.append('{0} {1}'.format(x,y))

  WktString += ','.join(Data)
  WktString += '))'

  return WktString

#-----------------------------------------------------------------------------------------

class Polygon():

  def __init__(self):

    self.x = []
    self.y = []

  #---------------------------------------------------------------------------------------

  def Load(self, XY):
    """Input polygon can be defined in two possible ways:
	
    1) list of x-y float pairs, e.g.
        [[22.0, -15.0],[24.0, -15.0],[24.0, -10.0],[22.0, -15.0]]
		
    2) wkt formatted string, e.g.
        'POLYGON((22. -15.,24. -15.,24. -10.,22. -15.))'

    """

    if type(XY) == list:
      # List of coordinate pairs
      self.x = [N[0] for N in XY]
      self.y = [N[1] for N in XY]

    elif type(XY) == str:
      # WKT String
      self.x, self.y = WktToXY(XY)

    else:
      print('Format not recognized')

  #---------------------------------------------------------------------------------------

  def Unwrap(self, Dir='Plus'):
    """
    """
    if Dir == 'Plus':
      self.x = [x if x > 0. else x+360. for x in self.x]
    else:
      self.x = [x if x < 0. else x-360. for x in self.x]

  #---------------------------------------------------------------------------------------

  def IsInside(self, x, y):

    x0 = self.x[0]
    y0 = self.y[0]

    inside = False
    n = len(self.x)

    for i in range(n+1):
      x1 = self.x[i % n]
      y1 = self.y[i % n]

      if min(y0,y1) < y <= max(y0,y1):
        if x <= max(x0,x1):
          if y0 != y1:
            xints = (y-y0)*(x1-x0)/(y1-y0)+x0
          if x0 == x1 or x <= xints:
            inside = not inside
      x0,y0 = x1,y1

    return inside

  #---------------------------------------------------------------------------------------

  def AddBuffer(self, Delta):
    """
    TO DO:
    Tempory implementation using Shapely.
    In the future, all Polygon objects will be defined this way

    """

    P = geometry.Polygon(zip(self.x, self.y))
    B = P.buffer(Delta)

    x, y = B.exterior.xy
    self.x = [i for i in x[:-1]]
    self.y = [i for i in y[:-1]]

  #---------------------------------------------------------------------------------------

  def Import (self, FileName, Type='xy'):

    with open(FileName, 'r') as f:

      XY = []
      if Type == 'wkt':
        XY = f.readline().strip()

      if Type == 'xy':
        for xy in f:
          xy = xy.strip()
          if xy:
            xy = re.split(',|;| ',xy)
            XY.append([float(xy[0]), float(xy[1])])

      self.Load(XY)
      f.close()
      return

    # Warn user if model file does not exist
    print('File not found.')

  #---------------------------------------------------------------------------------------

  def Area (self, Wgs=True):
    """
    Using Shoelace formula to compute area.
    Optionally, Wgs coordinates can be approximated to Km using
    sinusoidal projection (default).

    """

    if Wgs:
      x,y = WgsToXY(self.y, self.x)
    else:
      x = np.array(Lon)
      y = np.array(Lat)

    # Computing area
    A = np.dot(x, np.roll(y, 1))
    B = np.dot(y, np.roll(x, 1))

    return 0.5*np.abs(A-B)

  #---------------------------------------------------------------------------------------

  def CartGrid(self, Dx=0.1, Dy=0.1, Bounds=[]):
    """
    Produce a lat/lon cartesian grid.
    Dx and Dy distances are degrees (area is not preserved).
    Bounds are [MinX, MinY, MaxX, MaxY]

    """

    if Bounds:
      MinX = Bounds[0]
      MinY = Bounds[1]
      MaxX = Bounds[2]
      MaxY = Bounds[3]
    else:
      MinX = np.min(self.x)
      MinY = np.min(self.y)
      MaxX = np.max(self.x)
      MaxY = np.max(self.y)

    X = np.arange(MinX, MaxX, Dx)
    Y = np.arange(MinY, MaxY, Dy)

    XY = []
    for x in X:
      for y in Y:
        if self.IsInside(x, y):
          XY.append([x,y])

    return XY

  #---------------------------------------------------------------------------------------

  def SphereGrid(self, Delta=0.5, Bounds=[], Unwrap=False):
    """
    Distance between nearby points (in degree) is an approximated value.
    Bounds are [MinX, MinY, MaxX, MaxY]

    """

    X, Y = SphericalMesh(Delta)

    if Unwrap:
      X = [x+360. if x < 0. else x for x in X]

    if Bounds:
      MinX = Bounds[0]
      MinY = Bounds[1]
      MaxX = Bounds[2]
      MaxY = Bounds[3]
    else:
      MinX = np.min(self.x)
      MinY = np.min(self.y)
      MaxX = np.max(self.x)
      MaxY = np.max(self.y)

    XY = []
    for x, y in zip(X, Y):
      # Rough selection
      if x >= MinX and x <= MaxX:
        if y >= MinY and y <= MaxY:
          # Refined Selection
          if self.IsInside(x, y):
            XY.append([x,y])

    return XY

#-----------------------------------------------------------------------------------------

class Trace():

  def __init__(self):

    self.x = []
    self.y = []

  #---------------------------------------------------------------------------------------

  def Load(self, XY):
    """
    Input trace line can be defined in two possible ways:
	
    1) list of x-y float pairs, e.g.
        [[22.0, -15.0],[24.0, -15.0],[24.0, -10.0],[22.0, -15.0]]
		
    2) wkt formatted string, e.g.
        'LINESTRING((22. -15.,24. -15.,24. -10.,22. -15.))'

	"""

    if type(XY) == list:
      # List of coordinate pairs
      self.x = [N[0] for N in XY]
      self.y = [N[1] for N in XY]

    elif type(XY) == str:
      # WKT String
      self.x, self.y = WktToXY(XY)

    else:
      print('Format not recognized')

  #---------------------------------------------------------------------------------------

  def Buffer(self, Delta):
    """
    Return a polygon object containing the buffer area

    """

    L = geometry.LineString(zip(self.x, self.y))
    B = L.buffer(Delta)

    P = Polygon()
    x, y = B.exterior.xy

    P.x = [i for i in x[:-1]]
    P.y = [i for i in y[:-1]]

    return P

  #---------------------------------------------------------------------------------------

  def Resample(self, Delta):
    """
    Original code by Christian K (modified)
    https://stackoverflow.com/users/2588210/christian-k
    """

    x = np.array(self.x)
    y = np.array(self.y)

    xd = np.diff(x)
    yd = np.diff(y)

    dist = np.sqrt(xd**2+yd**2)

    u = np.cumsum(dist)
    u = np.hstack([[0],u])

    t = np.linspace(0,u.max(), int(u.max()//Delta))

    self.x = np.interp(t, u, x)
    self.y = np.interp(t, u, y)
