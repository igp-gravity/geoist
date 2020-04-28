# -*- coding: utf-8 -*-
"""
Created on Sun Apr 26 17:45:42 2020

@author: chens
"""

from geoist.inversion.regularization import Damping
from geoist.inversion.regularization import (Smoothness, Smoothness1D,
Smoothness2D, TotalVariation, TotalVariation1D, TotalVariation2D)
import numpy as np
from geoist.inversion.pfmodel import SmoothOperator

print('starting inversion')

shape = (2, 2)
regul = Damping(3) #最小构造约束
p = np.array([2, 3, 5])
print(regul.hessian(p).todense())
print(regul.gradient(p))

fd = np.array([[1, -1, 0],[0, -1, 1]]) #构造fd矩阵，包括深度加权在里面
smooth = Smoothness(fd)
smooth.value(p)
smooth.gradient(p)
smooth.hessian(p)

smooth1 = Smoothness1D(3)  #Smoothness的1维扩展
smooth1.hessian(p).todense()

smooth2 = Smoothness2D(shape) #Smoothness的2维扩展
smooth2.hessian(p).todense()

p1 = np.array([0, 0, 0])
tv1 = TotalVariation1D(2.67,3)
tv1.hessian(p).todense()

smop = SmoothOperator()
nz = 3
ny = 3
nx = 3
p = np.eye(nz*ny*nx).reshape(-1,nz,ny,nx)
sx = smop.derivation(p, component = 'dx').reshape(nz*ny*nx,-1)
sy = smop.derivation(p, component = 'dy').reshape(nz*ny*nx,-1)
sz = smop.derivation(p, component = 'dz').reshape(nz*ny*nx,-1)

ax = 0.1
ay = 0.2
az = 3.0

sm = np.vstack((ax*sx.T,ay*sy.T,az*sz.T))
smoothxyz = Smoothness(sm)
s1 = smoothxyz.hessian(np.ones(nx*ny*nz))
s11 = np.matmul(sm.T,sm)*2 #np.dot(sm.T,sm)*2

assert (s1==s11).all()

from geoist.inversion.hyper_param import LCurve
from geoist.pfm import inv3d
datamisfit = inv3d.Density3D(np.array(field.T).ravel(), [xp, yp, zp], mesh)
regul = smoothxyz.hessian(np.ones(nx*ny*nz)) #Damping(datamisfit.nparams)
regul_params = [10**i for i in range(-10, 5, 1)]
density3d = LCurve(datamisfit, regul, regul_params)
_ = density3d.fit()
print(density3d.regul_param_)
density3d.plot_lcurve()

sxtsx = smop.rderivation(smop.derivation(p, component = 'dx')).reshape(nz*ny*nx,-1)

sxx = smop.derivation(p, component = 'dxx').reshape(nz*ny*nx,-1)
syy = smop.derivation(p, component = 'dyy').reshape(nz*ny*nx,-1)
szz = smop.derivation(p, component = 'dzz').reshape(nz*ny*nx,-1)
sm2 = np.vstack((sxx.T,syy.T,szz.T))


