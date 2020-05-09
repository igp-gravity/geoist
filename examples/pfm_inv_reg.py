# -*- coding: utf-8 -*-
"""
Created on Sun Apr 26 17:45:42 2020

@author: chens
"""
import matplotlib.pyplot as plt
import numpy as np
from geoist import gridder
from geoist.inversion import geometry
from geoist.pfm import prism
from geoist.pfm import giutils
from geoist.inversion.mesh import PrismMesh
from geoist.vis import giplt
from geoist.inversion.regularization import Smoothness,Damping,TotalVariation
from geoist.inversion.pfmodel import SmoothOperator
from geoist.inversion.hyper_param import LCurve
from geoist.pfm import inv3d

meshfile = r"d:\msh.txt"
densfile = r"d:\den.txt"
#生成场源网格 NS-40km, EW-80km, 500个单元，z方向10层
area = (-20000, 20000, -40000, 40000, 2000, 32000) #NS EW Down
shape = (10, 20, 5) #nz ny nx
mesh = PrismMesh(area, shape)
density=np.zeros(shape)
density[3:8,9:12,1:4]=1.0 # z x y
mesh.addprop('density', density.ravel())
mesh.dump(meshfile, densfile, 'density') #输出网格到磁盘，MeshTools3D可视化
#生成核矩阵
kernel=[] 
narea = (-28750, 28750,-48750, 48750) #NS, EW
nshape = (30, 40) #NS, EW
depthz = []
xp, yp, zp = gridder.regular(narea, nshape, z=-1)
for i, layer in enumerate(mesh.layers()):
    for j, p in enumerate(layer):
        x1 = mesh.get_layer(i)[j].x1
        x2 = mesh.get_layer(i)[j].x2
        y1 = mesh.get_layer(i)[j].y1
        y2 = mesh.get_layer(i)[j].y2
        z1 = mesh.get_layer(i)[j].z1
        z2 = mesh.get_layer(i)[j].z2
        den = mesh.get_layer(i)[j].props
        model=[geometry.Prism(x1, x2, y1, y2, z1, z2, {'density': 1000.})]
        field = prism.gz(xp, yp, zp, model)
        kernel.append(field)       
        depthz.append((z1+z2)/2.0)
kk=np.array(kernel)        
kk=np.transpose(kernel)  #kernel matrix for inversion, 500 cells * 400 points
field0= np.mat(kk)*np.transpose(np.mat(density.ravel()))

field = field0 #giutils.contaminate(np.array(field0).ravel(), 0.05, percent = True)

#零均值
print(field.mean())
#field = field + 300. # field.mean()
#画图

plt.figure(figsize=(16, 8))
plt.subplot(1, 2, 1)
plt.title('gravity anomlay')
plt.axis('scaled')
levels = giplt.contourf(yp * 0.001, xp * 0.001, field0, nshape, 15)
cb = plt.colorbar(orientation='horizontal')
giplt.contour(yp * 0.001, xp * 0.001, field0, nshape,
                levels, clabel=False, linewidth=0.1)
plt.subplot(1, 2, 2)
plt.title('gravity anomlay with 5% noise')
plt.axis('scaled')
levels = giplt.contourf(yp * 0.001, xp * 0.001, field, nshape, 15)
cb = plt.colorbar(orientation='horizontal')
giplt.contour(yp * 0.001, xp * 0.001, field, nshape,
                levels, clabel=False, linewidth=0.1)
plt.show()

print('starting inversion')

smop = SmoothOperator()
nz = shape[0]
ny = shape[1]
nx = shape[2]
p = np.eye(nz*ny*nx).reshape(-1,nz,ny,nx)
sx = smop.derivation(p, component = 'dx').reshape(nz*ny*nx,-1)
sy = smop.derivation(p, component = 'dy').reshape(nz*ny*nx,-1)
sz = smop.derivation(p, component = 'dz').reshape(nz*ny*nx,-1)

am = 1000.0
ax = 300.0
ay = 300.0
az = 300.0
z0 = 10000 
beta = 1.0 #1.0

regul0 = Damping(nz*ny*nx)

wdepth = np.diag(1./(np.array(depthz)+z0)**beta)
sm0 = am*np.eye(nz*ny*nx)*wdepth
regul1 = Smoothness(sm0)

sm = np.vstack((am*np.eye(nz*ny*nx)*wdepth,
                az*np.dot(sz.T,wdepth),
                ay*np.dot(sy.T,wdepth),
                ax*np.dot(sx.T,wdepth)))
regul2 = Smoothness(sm) #np.eye(nz*ny*nx)

sm1 = np.vstack((az*np.dot(sz.T,wdepth),
                ay*np.dot(sy.T,wdepth),
                ax*np.dot(sx.T,wdepth)))
regul3 = Smoothness(sm1)

# L1 norm non-linear
regul4 = TotalVariation(0.0001, sm1)

#regul = Damping(nz*ny*nx)
datamisfit = inv3d.Density3D(np.array(field.T).ravel(), [xp, yp, zp], mesh
                             , movemean = False)

# 观测数据权
#weights = np.ones(datamisfit.nparams)
#datamisfit.set_weights(weights)

regul_params = [10**i for i in range(-5, 2, 1)]
# 目标函数
density3d = LCurve(datamisfit, regul4, regul_params, loglog=False)

initial = np.zeros(nz*ny*nx) #np.ones(datamisfit.nparams)
minval = initial 
maxval = initial + 1.0
bounds = list(zip(minval, maxval))
x0 = initial + 1.0
#_ = density3d.config('tcbound', bounds = bounds, nparams = datamisfit.nparams, x0  = x0).fit()
#_ = density3d.fit()
#solver.config('levmarq', initial=initial).fit()
_ = density3d.config('levmarq', initial=np.ones(datamisfit.nparams), maxit=100, maxsteps=100, tol=10**-4).fit()
print('Hyperparameter Lambda value is {}'.format(density3d.regul_param_))
density3d.plot_lcurve()

predicted = density3d[0].predicted()
residuals = density3d[0].residuals()

# print('bounds searching')
# solver = datamisfit + density3d.regul_param_*regul
# initial = np.zeros_like(density3d.p_) #np.ones(datamisfit.nparams)
# minval = initial 
# maxval = initial + 1.0
# bounds = np.zeros(2*len(initial))
# bounds[::2]=minval #奇数位置
# bounds[1::2]=maxval #偶数位置
#bounds = list(zip(minval, maxval))

#solver.config('acor', bounds = bounds, seed=0).fit()
#solver.config('acor', bounds = [0.5, 1.5], seed=0).fit()
#x0 = initial + 1.0
#solver.config('tcbound', bounds = bounds, x0  = x0).fit()
#solver.config('levmarq', initial=x0).fit()

# densinv = r"d:\deninv_solver.txt"
# values = np.fromiter(solver.estimate_, dtype=np.float)
# reordered = np.ravel(np.reshape(values, mesh.shape), order='F')
# np.savetxt(densinv, reordered, fmt='%.8f')    


plt.figure(figsize=(16, 8))
plt.subplot(1, 2, 1)
plt.axis('scaled')
plt.title('inversed gravity anomlay (mGal)')
levels = giplt.contourf(yp * 0.001, xp * 0.001, predicted, nshape, 15)
cb = plt.colorbar(orientation='horizontal')
giplt.contour(yp * 0.001, xp * 0.001, predicted, nshape,
                levels, clabel=False, linewidth=0.1)

plt.subplot(1, 2, 2)
plt.axis('scaled')
plt.title('residual (mGal)')
levels = giplt.contourf(yp * 0.001, xp * 0.001, residuals, nshape, 15)
cb = plt.colorbar(orientation='horizontal')
giplt.contour(yp * 0.001, xp * 0.001, residuals, nshape,
                levels, clabel=False, linewidth=0.1)

print('res mean={:.4f}; std={:.4f}'.format(residuals.mean(), residuals.std()))

densinv = r"d:\deninv_reg4.txt"
#values = np.fromiter(density3d.objectives[8].p_, dtype=np.float)
values = np.fromiter(density3d.p_, dtype=np.float)
reordered = np.ravel(np.reshape(values, mesh.shape), order='F')
np.savetxt(densinv, reordered, fmt='%.8f')    


