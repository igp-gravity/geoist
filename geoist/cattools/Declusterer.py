#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

import matplotlib.pyplot as plt
import numpy as np
import math as ma

import .Catalogue as Cat
import .CatUtils as CU

#-----------------------------------------------------------------------------------------

def GardnerKnopoff(M):

  Swin = 10**(0.1238*M+0.983)
  if M >= 6.5:
    Twin = 10.**(0.032*M+2.7389)*84600.
  else:
    Twin = 10.**(0.5409*M-0.547)*84600.
  return Twin, Swin

def Grunthal(M):

  Swin = np.exp(1.77+np.sqrt(0.037+1.02*M))
  if M >= 6.5:
    Twin = 10**(2.8+0.024*M)*84600.
  else:
    Twin = np.abs(np.exp(-3.95+np.sqrt(0.62+17.32*M)))*84600.
  return Twin, Swin

def Uhrhammer(M):

  Swin = np.exp(-1.024+0.804*M)
  Twin = np.exp(-2.87+1.235*M)*84600.
  return Twin, Swin

#-----------------------------------------------------------------------------------------

def WindowSearch(Db, WinFun=GardnerKnopoff, WinScale=1):
  """
  """

  def GetDate(Event):

    L = Event['Location'][0]
    S = CU.DateToSec(L['Year'],
                     L['Month'],
                     L['Day'],
                     L['Hour'],
                     L['Minute'],
                     L['Second'])
    return S

  def GetCoor(Event):
    X = Event['Location'][0]['Longitude']
    Y = Event['Location'][0]['Latitude']
    return [X, Y]

  def GetSize(Event):
    M = Event['Magnitude'][0]['MagSize']
    return M

  def DeltaSec(S0, S1):
    Sec = ma.fabs(S1-S0)
    return Sec

  def DeltaLen(C0, C1):
    Dis = CU.WgsDistance(C0[1],C0[0],C1[1],C1[0])
    return Dis

  def LogInfo(Type, Event, dT, dC):
    L = []
    L.append(Type)
    L.append(Event['Id'])
    L.append(Event['Location'][0]['Latitude'])
    L.append(Event['Location'][0]['Longitude'])
    L.append(Event['Magnitude'][0]['MagSize'])
    L.append(dT)
    L.append(dC)
    return L

  #---------------------------------------------------------------------------------------

  Events = []
  Log = []

  DbM = Db.Sort(Key='Magnitude', Owrite=0)
  Enum = DbM.Size()

  T0 = [0]*Enum
  S0 = [0]*Enum
  Tw = [0]*Enum
  Sw = [0]*Enum
  In = [0]*Enum

  for I in range(0,Enum):

    E0 = DbM.Events[I]
    T0[I] = GetDate(E0)
    S0[I] = GetCoor(E0)
    Tw[I], Sw[I] = WinFun(GetSize(E0))

  for I in range(0,Enum):
    if not In[I]:

      E0 = DbM.Events[I]
      Events.append(E0)
      Log.append([])
      Log[-1].append(LogInfo('MS', E0, 0, 0))

      for J in range(I+1,Enum):
        if not In[J]:

          E1 = DbM.Events[J]

          dT = DeltaSec(T0[I], T0[J])

          isAS = (T0[J] > T0[I] and dT <= Tw[I])
          isFS = (T0[J] < T0[I] and dT <= WinScale*Tw[I])

          if isAS or isFS:

            dS = DeltaLen(S0[I], S0[J])
            if dS <= Sw[I]:

              Flag = 'AS' if isAS else 'FS'
              Log[-1].append(LogInfo(Flag, E1, dT, dS))
              In[J] = 1

  DbM.Events = Events
  DbM.Sort()

  return DbM, Log

#-----------------------------------------------------------------------------------------

def PlotLog(Log, OutFile=[]):
  """
  Presently only for Debug
  """

  FS = [[],[]]
  AS = [[],[]]
  MS = [[],[]]

  fig = plt.figure(figsize=(5, 5))

  for L in Log:
    if len(L) > 1:
      for S in L:
        if S[0] == 'FS':
          FS[0].append(S[3])
          FS[1].append(S[2])
        if S[0] == 'AS':
          AS[0].append(S[3])
          AS[1].append(S[2])
        if S[0] == 'MS':
          MS[0].append(S[3])
          MS[1].append(S[2])

  plt.plot(AS[0], AS[1], 'b.', markersize=2)
  plt.plot(FS[0], FS[1], 'g.', markersize=2)
  plt.plot(MS[0], MS[1], 'r*', markersize=3)

  plt.xlabel('Longitude', fontsize=12, fontweight='bold')
  plt.ylabel('Latitude', fontsize=12, fontweight='bold')

  plt.grid('on')

  plt.tight_layout()

  plt.show(block=False)

  if OutFile:
    plt.savefig(OutFile, bbox_inches = 'tight', dpi = 150)