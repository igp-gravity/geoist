# -*- coding: utf-8 -*-
"""
Created on Fri Jan  4 16:55:10 2019
@author: chens
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
mesh = PrismMesh((0, 10000, 0, 20000, 0, 5000), (20, 20, 20))
mesh.addprop('density', 1000.0*np.random.rand(8000))
mesh.dump(meshfile, densfile, 'density') #输出网格到磁盘，MeshTools3D可视化
#生成核矩阵
kernel=[] 
xp, yp, zp = gridder.regular((-5000, 15000, -5000, 25000), (20, 20), z=-1)
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
kk=np.transpose(kk)  #kernel matrix for inversion, 500 cells * 400 points

#TO-DO list for Bei: understanding the forward and inversion of potential field
# 1. COO/CSR format for sparse matrix

# 2. SAVE the sensitivity matrix with HDF5 format

# 3. TEST compression by Haar wavelet approach

# 4. TEST compressed ratio and results

# 5. Comparison for JIT /numbaPRO/pypy or other potential CPU/GPU accelerate method

# 6. Add related codes to GEOIST library