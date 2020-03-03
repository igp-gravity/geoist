#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

import numpy as np
import scipy.odr as odr
import matplotlib.pyplot as plt

#-----------------------------------------------------------------------------------------

def PolyFun(B, X):
  """
  Arbitrary polynomial of arbitrary degree.
  Degree is taked from the length of the
  coefficient array (B).
  """

  X = np.array(X)
  Y = 0.0

  for d in np.arange(len(B)):
    Y += B[d]*(X**d)

  return Y

#-----------------------------------------------------------------------------------------

class Data(object):

  def __init__(self, X, Y, Xw=[], Yw=[]):

    self.X = np.array(X)
    self.Y = np.array(Y)

    self.Xw = np.array(Xw)
    self.Yw = np.array(Yw)

    if self.Xw.size == 0:
      self.Xw = np.zeros(self.X.size)
    if self.Yw.size == 0:
      self.Yw = np.zeros(self.Y.size)

    self.B = np.array([])
    self.E = np.array([])

  #-----------------------------------------------------------------------------------------

  def OrthReg(self, Deg=1):
    """
    Orthogonal regressiong of a polynomial of
    arbitrary degree.
    """

    B0 = np.ones(Deg+1)

    PF = odr.Model(PolyFun)

    if (self.Xw.all() == 0) and (self.Yw.all() == 0):
      Data = odr.Data(self.X, self.Y)
    else:
      Wx = 1./np.power(self.Xw,2)
      Wy = 1./np.power(self.Yw,2)
      Data = odr.Data(self.X, self.Y, Wx, Wy)

    Out = odr.ODR(Data, PF, beta0=B0)
    Out.run()

    self.B = Out.output.beta
    self.E = Out.output.sd_beta

  #-----------------------------------------------------------------------------------------

  def Plot(self, Xlabel='', Ylabel='', Axis=[], OutFile='', Dpi=300):
    """
    Plot utility to check regression results.
    """

    if self.X.size != 0 and self.Y.size != 0:

      plt.figure(figsize=(5,5))

      plt.errorbar(self.X, self.Y,
                   xerr=self.Xw, yerr=self.Yw,
                   fmt='o', capsize=5, elinewidth=1.0,
                   color=(0.4,0.4,0.4), ecolor=(0.8,0.8,0.8),
                   barsabove=True, label='Data')

      if not Axis:
        Axis = [np.min(self.X)-1,
                np.max(self.X)+1,
                np.min(self.Y)-1,
                np.max(self.Y)+1]

      Ax = np.linspace(Axis[0],Axis[1],100)

      plt.plot(Ax, Ax, color=(1,0,0),
                       linewidth=2,
                       label='1:1')

      if self.B.size != 0:
        plt.plot(Ax, PolyFun(self.B,Ax),
                 color=(0,0,0),
                 linewidth=2,
                 label='Regression')

      plt.legend(loc='upper left', numpoints=1)
      plt.xlabel(Xlabel)
      plt.ylabel(Ylabel)

      plt.grid(True)
      plt.axis(Axis)

      plt.show(block=False)

      if OutFile:
        plt.savefig(OutFile, bbox_inches='tight', dpi=Dpi)

  #-----------------------------------------------------------------------------------------

  def Close(self):

    plt.close('all')

  #-----------------------------------------------------------------------------------------

  def Print(self):
    """
    Print regression coefficients.
    """

    print('Coefficients: ['),
    print(', '.join([str(b) for b in self.B])),
    print(']\n')
