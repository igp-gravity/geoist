#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

import math as mt

from . import Regressor as Reg

#-----------------------------------------------------------------------------------------
# Generic 1:1

def Mw_Mw_Generic(MagSize, MagError):
  """
  Generic 1:1 conversion
  """

  return (MagSize, MagError)

#-----------------------------------------------------------------------------------------
# Generic from coefficients (polynomial)

def Mw_PolyReg(MagSize, MagError, Coeff):
  """
  Conversion using polynomial coefficients
  """

  M = Reg.PolyFun(Coeff, MagSize)
  E = None # Still to implement

  return (M, E)

#-----------------------------------------------------------------------------------------
# Scordilis (2006)
#
# MS >> Mw

def Ms_Mw_Scordilis2006(MagSize, MagError):
  """
  Linear
  """

  if MagSize >= 3.0 and MagSize <= 6.1:
    M = 0.67 * MagSize + 2.07
    E = mt.sqrt(0.17**2. + MagError**2.)

  elif MagSize > 6.1 and MagSize <= 8.2:
    M = 0.99 * MagSize + 0.08
    E = mt.sqrt(0.2**2. + MagError**2.)

  else:
    M = None
    E = None

  return (M, E)

# mb >> Mw

def mb_Mw_Scordilis2006(MagSize, MagError):
  """
  Linear
  """

  if MagSize >= 3.5 and MagSize <= 6.2:
    M = 0.85 * MagSize + 1.03
    E = mt.sqrt(0.29**2. + MagError**2.)

  else:
    M = None
    E = None

  return (M, E)

#-----------------------------------------------------------------------------------------
# Di Giacomo 2015
#
# Note:
#   No Error is provided from the paper.
#   Magnitude range is taken approximatively 
#   from the paper's pictures. To check...
#
# MS >> Mw

def Ms_Mw_Lin_DiGiacomo2015(MagSize, MagError):
  """
  Piecewise linear
  """

  if MagSize >= 3.5 and MagSize <= 6.47:
    M = 0.67 * MagSize + 2.13
    E = MagError

  elif MagSize > 6.47 and MagSize <= 8.0:
    M = 1.10 * MagSize - 0.67
    E = MagError

  else:
    M = None
    E = None

  return (M, E)

def Ms_Mw_Exp_DiGiacomo2015(MagSize, MagError):
  """
  Exponential
  """

  if MagSize >= 3.5 and MagSize <= 8.0:
    M = mt.exp(-0.222 + 0.233 * MagSize) + 2.863
    E = MagError

  else:
    M = None
    E = None

  return (M, E)

# mb >> Mw

def mb_Mw_Lin_DiGiacomo2015(MagSize, MagError):
  """
  Linear
  """

  if MagSize >= 4.0 and MagSize <= 6.5:
    M = 1.38 * MagSize - 1.79
    E = MagError

  else:
    M = None
    E = None

  return (M, E)

def mb_Mw_Exp_DiGiacomo2015(MagSize, MagError):
  """
  Exponential
  """

  if MagSize >= 4.0 and MagSize <= 6.5:
    M =  mt.exp(-4.664 + 0.859 * MagSize) + 4.555
    E = MagError

  else:
    M = None
    E = None

  return (M, E)

#-----------------------------------------------------------------------------------------
# Weatherill 2016
#
# MS >> Mw

def Ms_Mw_ISC_Weatherill2016(MagSize, MagError):
  """
  Piecewise linear
  """

  if MagSize >= 3.5 and MagSize <= 6.0:
    M = 0.616 * MagSize + 2.369
    E = mt.sqrt(0.147**2. + MagError**2.)

  elif MagSize > 6.0 and MagSize <= 8.0:
    M = 0.994 * MagSize + 0.1
    E = mt.sqrt(0.174**2. + MagError**2.)

  else:
    M = None
    E = None

  return (M, E)

def Ms_Mw_NEIC_Weatherill2016(MagSize, MagError):
  """
  Piecewise linear
  """

  if MagSize >= 3.5 and MagSize <= 6.47:
    M = 0.723 * MagSize + 1.798
    E = mt.sqrt(0.159**2. + MagError**2.)

  elif MagSize > 6.47 and MagSize <= 8.0:
    M = 1.005 * MagSize - 0.026
    E = mt.sqrt(0.187**2. + MagError**2.)

  else:
    M = None
    E = None

  return (M, E)

def Msz_Mw_NEIC_Weatherill2016(MagSize, MagError):
  """
  Piecewise linear
  """

  if MagSize >= 3.5 and MagSize <= 6.47:
    M = 0.707 * MagSize + 1.933
    E = mt.sqrt(0.179**2. + MagError**2.)

  elif MagSize > 6.47 and MagSize <= 8.0:
    M = 0.950 * MagSize + 0.359
    E = mt.sqrt(0.204**2. + MagError**2.)

  else:
    M = None
    E = None

  return (M, E)

# mb >> Mw

def mb_Mw_ISC_Weatherill2016(MagSize, MagError):
  """
  Linear
  """

  if MagSize >= 3.5 and MagSize <= 7.0:
    M = 1.048 * MagSize - 0.142
    E = mt.sqrt(0.317**2. + MagError**2.)

  else:
    M = None
    E = None

  return (M, E)

def mb_Mw_NEIC_Weatherill2016(MagSize, MagError):
  """
  Linear
  """

  if MagSize >= 3.5 and MagSize <= 7.0:
    M = 1.159 * MagSize - 0.659
    E = mt.sqrt(0.283**2. + MagError**2.)

  else:
    M = None
    E = None

  return (M, E)

#-----------------------------------------------------------------------------------------
# Edwards 2010
#
# Ml >> Mw

def Ml_Mw_Edwards2010(MagSize, MagError):
  """
  Polynomial
  """

  if MagSize <= 6.0:
    M = 0.0491 * MagSize**2 + 0.472 * MagSize + 1.02
    E = mt.sqrt(0.15**2. + MagError**2.)

  else:
    M = None
    E = None

  return (M, E)