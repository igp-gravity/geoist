# -*- coding: utf-8 -*-
"""
Created on Fri Jan  4 16:55:10 2019

@author: chens
"""
# 3rd imports
import matplotlib.pyplot as plt
import numpy as np
#from io import StringIO
# local imports
from geoist import gridder
from geoist.inversion import geometry
from geoist.pfm import prism
from geoist.inversion.mesh import PrismMesh
from geoist.vis import giplt
meshfile = r"d:\msh.txt" #StringIO()
densfile = r"d:\den.txt" #StringIO()
mesh = PrismMesh((0, 10, 0, 20, 0, 5), (5, 2, 2))
mesh.addprop('density', 1000.0*np.random.rand(20))
mesh.dump(meshfile, densfile, 'density')
#print(meshfile.getvalue().strip())
#print(densfile.getvalue().strip())
model=[]
for i, layer in enumerate(mesh.layers()):
    for j, p in enumerate(layer):
        #print(i,j, p)
        x1 = mesh.get_layer(i)[j].x1
        x2 = mesh.get_layer(i)[j].x2
        y1 = mesh.get_layer(i)[j].y1
        y2 = mesh.get_layer(i)[j].y2
        z1 = mesh.get_layer(i)[j].z1
        z2 = mesh.get_layer(i)[j].z2
        den = mesh.get_layer(i)[j].props
        model.append(geometry.Prism(x1, x2, y1, y2, z1, z2, den))
        #model.append(geometry.Prism(x1, x2, y1, y2, z1, z2, {'density': 1000}))

xp, yp, zp = gridder.regular((-5, 15, -5, 25), (20, 20), z=-1)
field = prism.gz(xp, yp, zp, model)

plt.figure(figsize=(8, 9))
plt.axis('scaled')
plt.title('forward gravity anomaly')
levels = giplt.contourf(yp , xp , field, (20, 20), 15)
cb = plt.colorbar()
giplt.contour(yp , xp , field, (20, 20),
                levels, clabel=False, linewidth=0.1)
plt.show()