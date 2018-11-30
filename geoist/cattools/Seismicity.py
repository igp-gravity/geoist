#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import scipy.optimize as spo
import matplotlib.pyplot as plt

import .Selection as Sel

#-----------------------------------------------------------------------------------------

def TableMerge(M0, M1, Y0, Y1):
  """
  Utility to assemble a completeness table from arrays
  """

  CompTable = []

  for m0, m1, y0, y1 in zip(M0, M1, Y0, Y1):
    CompTable.append([m0, m1, y0, y1])

  return CompTable


def TableSplit(CompTable):
  """
  Utility to split a completeness table into arrays
  """

  M0 = [m[0] for m in CompTable]
  M1 = [m[1] for m in CompTable]
  Y0 = [m[2] for m in CompTable]
  Y1 = [m[3] for m in CompTable]

  return M0, M1, Y0, Y1

#-----------------------------------------------------------------------------------------

def GetEventRates(Db, CompTable, Area=1.):
  """
  Method to compute observed annual rates (incremental and cumulative) from a given
  completeness table. In this implementation, completeness is one window per
  magnitude bin in M. Example:
  CompTable = [[4.50, 0.25, 2000., 2013.],
               [4.75, 0.25, 1980., 2013.],
               [5.00, 0.25, 1970., 2013.],
               [5.25, 0.25, 1960., 2013.],
               [5.50, 0.50, 1950., 2013.],
               [6.00, 1.50, 1901., 2013.]]
  """

  Enum = []
  Data = [[],[]]

  for CT in CompTable:

    MinM = CT[0]
    MaxM = CT[0]+CT[1]
    MinY = CT[2]
    MaxY = CT[3]

    # Catalogue selection (Magnitude-Year)
    DbM = Sel.MagRangeSelect(Db, MinM, MaxM)
    DbY = Sel.TimeSelect(DbM, MinY, MaxY)

    # Computing incremental rates
    RY = float(DbY.Size())/float(MaxY-MinY)
    Enum.append(RY/float(Area))

    # Data per magnitude bin
    Data[0].append(DbY.Extract(Key='Year'))
    Data[1].append(DbY.Extract(Key='MagSize'))

  # Cumulative distribution
  Ecum = np.cumsum(Enum[::-1])[::-1]

  return Enum, Ecum, Data

#-----------------------------------------------------------------------------------------

def MfdCum(a, b, Mbin, Mmax):
  """
  Cumulative MFD (Truncated Gutenberg-Richter)
  """

  Mbin = np.array(Mbin)

  Enum = (10.**a)*(10.**(-b*Mbin)-10.**(-b*Mmax))

  # If M > Mmax, remove negative rates
  Enum[Enum < 0] = 0

  return Enum

#-----------------------------------------------------------------------------------------

def MfdInc(a, b, Mbin, Minc, Mmax):
  """
  Incremental MFD (for discrete magnitude intervals).
  """

  Mbin = np.array(Mbin)
  Minc = np.array(Minc)

  Enum0 = MfdCum(a, b, Mbin, Mmax)
  Enum1 = MfdCum(a, b, Mbin+Minc, Mmax)

  return (Enum0-Enum1)

#-----------------------------------------------------------------------------------------

def MfdFit(ab, Enum, Mbin, Minc, Mmax, Merr, bfix=[]):
  """
  Misfit function (log normal)
  """

  # Target coefficients
  a = ab[0]

  if not bfix:
    b = ab[1]
  else:
    b = bfix

  # Analytical distribution
  Esyn = MfdInc(a, b, Mbin, Minc, Mmax)

  # Avoid log of 0
  Esyn[Esyn <= 0] = 1e-300

  Eres = np.log10(Enum/Esyn)

  # L2-norm
  Mfit = np.sum((Eres/Merr)**2.)

  return Mfit

#-----------------------------------------------------------------------------------------

def MfdOptimize(Enum, Mbin, Minc, Mmax, Merr=[], a0=[], b0=[], bfix=[]):
  """
  Optimisation function
  Note: Minc and Merr can be single (constant) values or array
  """

  # Convert to numpy array
  Enum = np.array(Enum)
  Mbin = np.array(Mbin)
  Minc = np.array(Minc)
  Merr = np.array(Merr)

  # Setting initial values for the search
  if not a0: a0 = 10.
  if not b0: b0 = 1.

  if Merr.size == 0:
    Merr = np.ones_like(Enum)
  if Merr.size == 1:
    Merr = Merr * np.ones_like(Enum)
  if Minc.size == 1:
    Minc = Minc * np.ones_like(Enum)

  # Remove zero rate bins
  idx = (Enum > 0.0)
  Enum = Enum[idx]
  Mbin = Mbin[idx]
  Minc = Minc[idx]
  Merr = Merr[idx]

  # Optimization of GR coefficients
  Out = spo.minimize(MfdFit, [a0, b0], args=(Enum, Mbin, Minc, Mmax, Merr, bfix))

  a = Out.x[0]
  b = Out.x[1]

  if bfix:
    b = bfix

  return a, b

#-------------------------------------------
# Plot results

def MfdPlot(a, b, Mmax, Enum=[], Ecum=[], Mbin=[], Minc=[], OutFile=[]):

  # Convert to numpy array
  Enum = np.array(Enum)
  Mbin = np.array(Mbin)
  Minc = np.array(Minc)

  # Plot
  plt.figure(figsize=(6,4))

  # Observed Incremental rates
  if any(Enum) and any(Mbin) and any(Minc):
    plt.bar(Mbin, Enum, Minc, edgecolor=[[0,0,0] for i in range(0,len(Mbin))],
                              color=[0.9,0.9,0.9],
                              linewidth=1,
                              label='Observed Incremental',
                              align='edge',
                              log=True,
                              zorder=1)

  # Observed Cumulative rates
  if any(Enum) and any(Mbin):
    plt.plot(Mbin, Ecum, 'o', color=[1,0,0],
                              markersize=6,
                              markeredgecolor=[0,0,0],
                              linewidth=2,
                              label='Observed Cumulative',
                              zorder=4)

  # Inverted Incremental rates
  Ninc = MfdInc(a, b, Mbin, Minc, Mmax)
  plt.plot(Mbin+Minc/2., Ninc, 's', color=[1,1,1],
                                    markersize=8,
                                    markeredgecolor=[0,0,0],
                                    linewidth=2,
                                    label='Inverted Incremental',
                                    zorder=3)

  # Inverted Cumulative rates
  Maxs = np.arange(min(Mbin), Mmax, 0.0001)
  Ncum = MfdCum(a, b, Maxs, Mmax)
  plt.plot(Maxs, Ncum, color=[1,0,0],
                       linewidth=2,
                       label='Inverted Cumulative',
                       zorder=2)

  plt.title('Truncated G-R Distribution')
  plt.legend(loc=1, borderaxespad=0., numpoints=1)

  plt.xlabel('Magnitude')
  plt.ylabel('Occurrence Rate (Event/Year)')

  plt.gca().yaxis.grid(color='0.65',linestyle='--')
  plt.tight_layout()

  plt.xlim((min(Mbin)-0.5, Mmax+0.5))
  plt.ylim((0.01*min(Enum), 10*max(Ncum)))

  plt.show(block=False)

  if OutFile:
    plt.savefig(OutFile, bbox_inches='tight', dpi=150)

