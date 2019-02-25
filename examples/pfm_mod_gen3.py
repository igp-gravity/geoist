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
from geoist.inversion.regularization import Damping

meshfile = r"d:\msh.txt"
densfile = r"d:\den.txt"
#生成场源网格 10km*20km*5km, 500个单元，z方向5层
area = (-47500, 47500, -47500, 47500, 4000, 32000)
shape = (20, 20, 5)
mesh = PrismMesh(area, shape)
density=np.zeros(shape)
density[9:12,9:12,2:3]=1.00
mesh.addprop('density', density.ravel())
mesh.dump(meshfile, densfile, 'density') #输出网格到磁盘，MeshTools3D可视化
#生成核矩阵
kernel=[] 
narea = (-48750, 48750,-48750, 48750)
nshape = (40, 40)
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
        model=[geometry.Prism(x1, x2, y1, y2, z1, z2, {'density': 1000})]
        field = prism.gz(xp, yp, zp, model)
        kernel.append(field)       
kk=np.array(kernel)        
kk=np.transpose(kernel)  #kernel matrix for inversion, 500 cells * 400 points
field=np.mat(kk)*np.transpose(np.mat(density.ravel()))
#画图
plt.title('gz')
levels = giplt.contourf(yp * 0.001, xp * 0.001, field, nshape, 15)
cb = plt.colorbar()
giplt.contour(yp * 0.001, xp * 0.001, field, shape,
                levels, clabel=False, linewidth=0.1)
plt.show()

###反演

#datamisfit = (1./tts.size)*srtomo.SRTomo(tts, srcs, recs, mesh)
#regul = Damping(mesh.shape)
#regul_params = [10**i for i in range(-10, -2, 1)]
#tomo = LCurve(datamisfit, regul, regul_params)
# _ = tomo.fit()
#print(numpy.array_repr(tomo.estimate_.reshape(shape), precision=0))
#
#predicted = tomo[0].predicted()
#residuals = tomo[0].residuals()
#print '%.4f %.4f' % (residuals.mean(), residuals.std())
#tomo.regul_param_
#par_tomo = LCurve(datamisfit, regul, regul_params, njobs=2)
#_ = par_tomo.fit()  # Will you 2 processes to run inversions
#par_tomo.regul_param_
# initial = numpy.ones(mesh.size)
#_ = tomo.config('newton', initial=initial, tol=0.2).fit()
#tomo.regul_param_








