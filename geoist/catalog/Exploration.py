#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as sig

from . import Selection as Sel

#-----------------------------------------------------------------------------------------

def GetHypocenter(Db, All=False):

  x = Db.Extract('Longitude', All)
  y = Db.Extract('Latitude', All)
  z = Db.Extract('Depth', All)

  return x, y, z

#-----------------------------------------------------------------------------------------

def GetMagnitudePair(Db, Code1, Code2):

  Mout = [[],[],[],[]]

  for E in Db.Events:
    m1 = None
    m2 = None

    for M in E['Magnitude']:
      MC = M['MagCode']
      MT = M['MagType']
      MS = M['MagSize']
      ME = M['MagError']

      if (MC == Code1[0] or Code1[0] == '*') and MT == Code1[1]:
        m1 = MS
        e1 = ME
      if (MC == Code2[0] or Code2[0] == '*') and MT == Code2[1]:
        m2 = MS
        e2 = ME

    if m1 and m2:
      Mout[0].append(m1)
      Mout[1].append(m2)
      Mout[2].append(e1)
      Mout[3].append(e2)

  return Mout

#-----------------------------------------------------------------------------------------

def GetKeyHisto(Db, Key, Bins=[], Bmin=[], Bmax=[], Bnum=10, Blog=False,
                         Norm=True, Plot=True, OutFile=[]):

  Data = Db.Extract(Key)

  # Remove Nans
  Data = [D for D in Data if D is not None]

  if not Bins:
    if not Bmin:
      Bmin = min(Data)
    if not Bmax:
      Bmax = max(Data)
    if Blog:
      Bins = np.logspace(np.log10(Bmin), np.log10(Bmax), Bnum)
    else:
      Bins = np.linspace(Bmin, Bmax, Bnum)

  Hist = np.histogram(Data, Bins)[0]
  Bmid = np.diff(Bins)/2.+Bins[:-1]
  Bdlt = np.diff(Bins)

  if Norm:
   Hist = Hist.astype('float32') / len(Data)

  # Plot time histogram
  if Plot:
    fig = plt.figure(figsize=(6,3.5))

    plt.bar(Bmid, Hist, Bdlt, color=[1,0,0], edgecolor=[1,1,1])

    plt.xlabel(Key, fontsize=14, fontweight='bold')
    plt.ylabel('Nr. Events', fontsize=14, fontweight='bold')

    plt.tight_layout()
    plt.show(block=False)

    if OutFile:
      plt.savefig(OutFile, bbox_inches='tight', dpi=150)

  return Hist, Bmid

#-----------------------------------------------------------------------------------------

def AgencyReport(Db, Code, Key=[], LogFile=[], Threshold=0):

  if Code in ['Magnitude','Mag','M']:
    ItL, ItD = Db.KeyStat('MagCode')
  elif Code in ['Location','Loc','L']:
    ItL, ItD = Db.KeyStat('LocCode')
  else:
    print('Error: No valid code')
    return

  # Report only specific keys
  if Key:
    ItLs = []
    ItDs = {}

    if type(Key) != list:
      Key = [Key]

    for K in Key:
     if K in ItL:
       ItLs.append(K)
       ItDs[K] = ItD[K]
    ItL = ItLs
    ItD = ItDs

  StrLog = ''

  for It in ItL:
    if ItD[It] >= Threshold:
      StrLog += 'Agency: {0} | Occurrence: {1}'.format(It, ItD[It])

      if Code in ['Magnitude','Mag','M']:
        DbC = Db.Filter('MagCode',It,Owrite=False)
        MaL, MaD = DbC.KeyStat('MagType')

        StrLog += ' | Types:'
        for Ma in MaL:
          StrLog += ' {0} ({1})'.format(Ma, MaD[Ma])

      StrLog += '\n'
    else:
      break

  if LogFile:
    # Open input ascii file
    with open(LogFile, 'w') as f:
      f.write(StrLog)
      f.close()
      return
    # Warn user if model file does not exist
    print('Cannot open file')

  else:
    print(StrLog)

#-----------------------------------------------------------------------------------------

def KeyTimeHisto(Db, Code, Key=[],
                           Year0=[], Year1=[], Delta=5,
                           Threshold=0, OutFile=[]):

  if not Year0:
    Year0 = min(Db.Extract('Year'))
  if not Year1:
    Year1 = max(Db.Extract('Year'))
  YBins = np.arange(Year0, Year1+Delta, Delta)

  ItL, ItD = Db.KeyStat(Code)

  # Filter by threshold
  ItL = [K for K in ItL if ItD[K] >= Threshold]
  ItD = {K:V for (K,V) in ItD.items() if V >= Threshold}

  # Filter by key
  if Key:
    ItL = [K for K in ItL if K in Key]
    ItD = {K:ItD[K] for K in ItL}

  for N, Agn in enumerate(ItL):

    DbA = Db.Filter(Code, Agn, Owrite=0)
    YearArray = DbA.Extract('Year')
    NewRow = np.histogram(YearArray, YBins)

    if N == 0:
      Histo = NewRow[0]
    else:
      Histo = np.vstack([Histo, NewRow[0]])

  # Plot time histogram
  fig = plt.figure(figsize=(8, 5))

  X = YBins
  Y = np.arange(0, len(ItL)+1)
  Z = np.log(Histo.clip(min=1E-10))

  plt.pcolor(X, Y, Z, cmap='Purples',
                      vmin=0,
                      vmax=np.max(Z))

  plt.xticks(X, map(str,X), rotation='45')
  plt.yticks(Y+0.5, ItL, rotation='horizontal')
  plt.margins(0)

  plt.gca().yaxis.tick_right()
  plt.axes().yaxis.grid(True)

  plt.gca().xaxis.set_ticks_position('none')
  plt.gca().yaxis.set_ticks_position('none')

  plt.xlabel('Year', fontsize=14, fontweight='bold')
  plt.ylabel('Agency Code', fontsize=14, fontweight='bold')

  plt.tight_layout()

  plt.show(block=False)

  if OutFile:
    plt.savefig(OutFile, bbox_inches='tight', dpi=150)

#-----------------------------------------------------------------------------------------

def MagTimeBars(Db, Mag0=[], Mag1=[], MBin=0.5,
                    Year0=[], Year1=[], Delta=5,
                    OutFile=[]):

  if not Mag0:
    Mag0 = min(Db.Extract('MagSize'))
  if not Mag1:
    Mag1 = max(Db.Extract('MagSize'))
  MBins = np.arange(Mag0, Mag1+MBin, MBin)

  if not Year0:
    Year0 = min(Db.Extract('Year'))
  if not Year1:
    Year1 = max(Db.Extract('Year'))
  YBins = np.arange(Year0, Year1+Delta, Delta)

  plt.figure(figsize=(8, 4))

  for C,MB in enumerate(MBins):

    DbM = Db.Filter('MagSize', MB, Opr='>=', Owrite=0)
    YArray = DbM.Extract('Year')
    YHist = np.histogram(YArray, YBins)[0]

    Cnum = float(len(MBins))
    C = (Cnum-C)/Cnum

    X = YBins[:-1]
    Y = YHist

    if any(Y):
      plt.bar(X, Y, Delta, color=[C,C,C],
                           log=True,
                           label=r'$\geq${0}'.format(MB))

  plt.xticks(X, map(str,X), rotation='45')
  plt.margins(0)

  plt.gca().yaxis.tick_right()

  plt.gca().xaxis.set_ticks_position('none')
  plt.gca().yaxis.set_ticks_position('none')

  plt.xlabel('Years', fontsize=14, fontweight='bold')
  plt.ylabel('Nr. Events', fontsize=14, fontweight='bold')

  plt.tight_layout()
  plt.legend(loc=2)
  plt.show(block=False)

  if OutFile:
    plt.savefig(OutFile, bbox_inches='tight', dpi=150)

#-----------------------------------------------------------------------------------------

def MagTimePlot(Db, Mag0=[], Mag1=[],
                    Year0=[], Year1=[],
                    CompTable=[], 
                    OutFile=[]):

  if not Mag0:
    Mag0 = min(Db.Extract('MagSize'))
  if not Mag1:
    Mag1 = max(Db.Extract('MagSize'))

  if not Year0:
    Year0 = min(Db.Extract('Year'))
  if not Year1:
    Year1 = max(Db.Extract('Year'))

  DbS = Sel.MagRangeSelect(Db, Mag0, Mag1, Owrite=0, TopEdge=True)
  DbS = Sel.TimeSelect(DbS, Year0, Year1, Owrite=0)

  X = DbS.Extract('Year')
  Y = DbS.Extract('MagSize')

  plt.figure(figsize=(7, 4))

  plt.plot(X, Y, 'o',markersize=3,
                  color=[0,0,0],
                  markeredgecolor=[0,0,0],
                  markeredgewidth=1.5)

  # Plot completeness
  if CompTable:
    PlotCompTable(CompTable)

  plt.gca().yaxis.grid(color='0.',linestyle='-')
  plt.gca().xaxis.grid(color='0.65',linestyle='--')

  plt.gca().xaxis.set_ticks_position('none')
  plt.gca().yaxis.set_ticks_position('none')

  plt.title('Time-Magnitude Distribution', fontsize=14, fontweight='bold')
  plt.xlabel('Years', fontsize=14, fontweight='bold')
  plt.ylabel('Magnitude', fontsize=14, fontweight='bold')

  plt.axis([Year0, Year1, Mag0, Mag1])
  # plt.tight_layout()
  plt.show(block=False)

  if OutFile:
    plt.savefig(OutFile, bbox_inches='tight', dpi=150)

#-----------------------------------------------------------------------------------------

def RateDensityPlot(Db, Mag0=[], Mag1=[], MBin=0.25,
                        Year0=[], Year1=[], Delta=2,
                        CompTable=[], 
                        Normalise=True,
                        OutFile=[]):

  if not Mag0:
    Mag0 = min(Db.Extract('MagSize'))
  if not Mag1:
    Mag1 = max(Db.Extract('MagSize'))
  MBins = np.arange(Mag0, Mag1+MBin, MBin)

  if not Year0:
    Year0 = min(Db.Extract('Year'))
  if not Year1:
    Year1 = max(Db.Extract('Year'))
  YBins = np.arange(Year0, Year1+Delta, Delta)

  Histo = np.zeros((np.size(MBins), np.size(YBins)))

  # Catalogue selection (Magnitude-Year)
  DbM = Sel.MagRangeSelect(Db, Mag0, Mag1, TopEdge=True)
  DbY = Sel.TimeSelect(DbM, Year0, Year1)

  M = DbY.Extract('MagSize')
  Y = DbY.Extract('Year')

  Hist = np.histogram2d(Y, M, bins=(YBins, MBins))[0]
  Hist = np.transpose(Hist)

  if Normalise:
    for I in range(0,len(Hist)):
      Max = np.max(Hist[I])
      if Max > 0:
        Hist[I] = Hist[I]/Max

  # Plot
  plt.figure(figsize=(7, 4))

  plt.pcolormesh(YBins, MBins, Hist, cmap='Greys', vmin=0)

  # Plot completeness
  if CompTable:
    PlotCompTable(CompTable)

  plt.gca().xaxis.grid(color='0.65',linestyle='--')
  plt.gca().yaxis.grid(color='0.',linestyle='-')

  plt.gca().xaxis.set_ticks_position('none')
  plt.gca().yaxis.set_ticks_position('none')

  plt.title('Occurrence Rate Density', fontsize=14, fontweight='bold')
  plt.xlabel('Years', fontsize=12, fontweight='bold')
  plt.ylabel('Magnitude', fontsize=12, fontweight='bold')

  plt.gca().xaxis.grid(color='0.65',linestyle='-')
  plt.gca().yaxis.grid(color='0.65',linestyle='-')

  plt.axis([Year0, Year1, Mag0, Mag1])
  # plt.tight_layout()
  plt.show(block=False)

  if OutFile:
    plt.savefig(OutFile, bbox_inches = 'tight', dpi = 150)


def PlotCompTable(CompTable):

  for CT in CompTable:

    X = [CT[2], CT[3], CT[3], CT[2], CT[2]]
    Y = [CT[0], CT[0], CT[0]+CT[1], CT[0]+CT[1], CT[0]]

    plt.plot(X, Y, 'r--', linewidth=2)
    plt.fill(X, Y, color='y',alpha=0.1)

#-----------------------------------------------------------------------------------------

def DuplicateCheck(Log, Tmax=[], Smax=[],
                        Tnum=[], Snum=[],
                        Smooth=[],
                        OutFile=[]):
  """
  """

  dT = [I[4] for I in Log if I[4] > 0]
  dS = [I[5] for I in Log if I[5] > 0]

  if not Tmax:
    Tmax = np.max(dT)
  if not Smax:
    Smax = np.max(dS)
  if not Tnum:
    Tnum = 100
  if not Snum:
    Snum = 100

  XBins = np.linspace(0, Tmax, Tnum)
  YBins = np.linspace(0, Smax, Snum)
  
  if(len(dT)-len(dS)>0):
      dT=dT[:len(dS)]
  elif(len(dT)-len(dS)<0):
      dS=dS[:len(dT)]
      
  H = np.histogram2d(dT, dS, [YBins, XBins])

  def Gaussian(Size,Sigma):
    x = np.arange(0, Size[0], 1, float)
    y = np.arange(0, Size[1], 1, float)
    Gx = np.exp(-(x-Size[0]/2)**2/Sigma[0]**2)
    Gy = np.exp(-(y-Size[1]/2)**2/Sigma[1]**2)
    return np.outer(Gy,Gx)

  if any(Smooth):
    # kern = np.ones((Smooth,Smooth))/Smooth
    kern = Gaussian((Tnum,Snum),Smooth)
    H0 = sig.convolve2d(H[0], kern, mode='same')
  else:
    H0 = H[0]

  # Plot time histogram
  fig = plt.figure(figsize=(5, 5))

  plt.pcolor(XBins, YBins, H0, cmap='Purples')
  plt.xlabel('Time', fontsize=12, fontweight='bold')
  plt.ylabel('Distance', fontsize=12, fontweight='bold')
  plt.grid('on')
  plt.tight_layout()
  plt.show(block=False)

  if OutFile:
    plt.savefig(OutFile, bbox_inches='tight', dpi=150)