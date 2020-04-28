# -*- coding: utf-8 -*-
"""
Created on Sun Jan  6 14:28:35 2019

@author: lhl
"""
import matplotlib.pyplot as plt
import numpy as np
from geoist import gridder
from geoist.inversion import geometry
from geoist.pfm import prism
from geoist.inversion.mesh import PrismMesh
from geoist.vis import giplt


meshfile = r"d:\msh.txt"
densfile = r"d:\den.txt"
#生成场源网格 10km*20km*5km, 500个单元，z方向5层
area = (-20000, 20000, -40000, 40000, 2000, 32000) #ny nx nz
shape = (10, 20, 5) #nz nx ny
mesh = PrismMesh(area, shape)
density=np.zeros(shape)
density[3:8,9:12,1:4]=1.0 # z x y
mesh.addprop('density', density.ravel())
mesh.dump(meshfile, densfile, 'density') #输出网格到磁盘，MeshTools3D可视化
#生成核矩阵
kernel=[] 
narea = (-28750, 28750,-48750, 48750)
nshape = (30, 40)
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
kk=np.array(kernel)        
kk=np.transpose(kernel)  #kernel matrix for inversion, 500 cells * 400 points
field=np.mat(kk)*np.transpose(np.mat(density.ravel()))
#画图
plt.figure()
plt.title('gz')
plt.axis('scaled')
levels = giplt.contourf(yp * 0.001, xp * 0.001, field, nshape, 15)
cb = plt.colorbar()
giplt.contour(yp * 0.001, xp * 0.001, field, nshape,
                levels, clabel=False, linewidth=0.1)
plt.show()

###反演
print('starting inversion')
from geoist.inversion.regularization import Damping
from geoist.inversion.hyper_param import LCurve
from geoist.pfm import inv3d
datamisfit = inv3d.Density3D(np.array(field.T).ravel(), [xp, yp, zp], mesh)
regul = Damping(datamisfit.nparams)
regul_params = [10**i for i in range(-10, 1, 1)]
density3d = LCurve(datamisfit, regul, regul_params, loglog=False)
_ = density3d.fit()
print(density3d.regul_param_)
density3d.plot_lcurve()
#print(np.array_repr(density3d.estimate_.reshape(shape), precision=0))
#density3d = LCurve(datamisfit, regul, regul_params, njobs=4)

#density3d = LCurve(datamisfit, regul, regul_params, njobs=2)
#_ = density3d.fit()  # Will you 2 processes to run inversions
#density3d.regul_param_
# initial = numpy.ones(mesh.size)
#_ = density3d.config('newton', initial=initial, tol=0.2).fit()
#density3d.regul_param_


predicted = density3d[0].predicted()
residuals = density3d[0].residuals()


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

densinv = r"d:\deninv1.txt"
values = np.fromiter(density3d.p_, dtype=np.float)
reordered = np.ravel(np.reshape(values, mesh.shape), order='F')
np.savetxt(densinv, reordered, fmt='%.8f')    

normfile = r"d:\norms.txt"
with open(normfile, 'w') as f:
    for i in range(len(density3d.dnorm)):
        f.write('{} {} {}\n'.format(density3d.mnorm[i], density3d.dnorm[i]
                                    ,regul_params[i]))








