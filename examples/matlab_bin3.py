# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 21:33:32 2020

@author: chens
"""
import scipy.io as sio
import matplotlib.pyplot as plt
import numpy as np
from geoist import gridder
from geoist.inversion import geometry
from geoist.pfm import prism, giutils
from geoist.inversion.mesh import PrismMesh
from geoist.vis import giplt

#fname1 = "D:\\projects\\model1-2-to-Jin-Chen\\rho1.mat"
fname2 = "D:\\projects\\model1-2-to-Jin-Chen\\rho3.mat"
data2=sio.loadmat(fname2)
#data2=sio.loadmat(fname2)
rho10 = data2['data_dens3']
#rho2 = data2['data_dens2']
rho1 = rho10[::4,::10]

maglist = [10.0, 6.0, 2.0, 6.0, 30., 20.0, 40., 4.0, 50.]
mag1 = rho1.copy()
for i, di in enumerate(set(rho1.ravel())):
    print(di, maglist[i])
    mag1[rho1 == di] = maglist[i]
    
    
plt.figure()
plt.imshow(rho1)

plt.figure()
plt.imshow(mag1)

#plt.figure()
#plt.imshow(rho2) #, origin='lower')

meshfile = r"d:\msh3.txt"
densfile = r"d:\den3.txt"
magfile = r"d:\mag3.txt"
graoutfile = r"d:\gra3.dat"
magoutfile = r"d:\mag3.dat"
graoutfile1 = r"d:\gra3n.dat"
magoutfile1 = r"d:\mag3n.dat"
area = (-200, 200, -600, 600, 0, 600) #x y z
shape = (100, 150, 1) # z y x
mesh = PrismMesh(area, shape)
mesh.addprop('density', 1000.*rho1.ravel()-2650.0001)
mesh.addprop('magnetization', mag1.ravel())
mesh.dump(meshfile, densfile, 'density') #输出网格到磁盘，MeshTools3D可视化
mesh.dump(meshfile, magfile, 'magnetization')
#生成核矩阵
kernel=[] 
narea = (-500, 500,-1000, 1000) #y x
nshape = (20, 40)
xp, yp, zp = gridder.regular(narea, nshape, z=-1)
prisms=[]
for p in mesh:
    prisms.append(p)
print('kernel')
inc, dec = 30, -4
kernelgz = prism.gz_kernel(xp, yp, zp, prisms)
for i, layer in enumerate(mesh.layers()):
    for j, p in enumerate(layer):
        x1 = mesh.get_layer(i)[j].x1
        x2 = mesh.get_layer(i)[j].x2
        y1 = mesh.get_layer(i)[j].y1
        y2 = mesh.get_layer(i)[j].y2
        z1 = mesh.get_layer(i)[j].z1
        z2 = mesh.get_layer(i)[j].z2
        den = mesh.get_layer(i)[j].props
        model=[geometry.Prism(x1, x2, y1, y2, z1, z2, 
                              {'magnetization': giutils.ang2vec(1, inc, dec)})]
        field = prism.tf(xp, yp, zp, model, inc, dec)
        kernel.append(field)       
     
kk=np.transpose(kernel)  #kernel matrix for inversion, 500 cells * 400 points
field_mag=np.mat(kk)*np.transpose(np.mat(mag1.ravel()))
field_gra=np.mat(kernelgz)*np.transpose(np.mat(rho1.ravel()))
field_mag1 = giutils.contaminate(np.array(field_mag).ravel(), 0.05, percent = True)
field_gra1 = giutils.contaminate(np.array(field_gra).ravel(), 0.05, percent = True)

#保存正演异常
with open(graoutfile, 'w') as f:
    f.write('! model 2 gravity anomlay (mGal)\n')
    f.write('{}\n'.format(len(field_gra)))
    for i in range(len(field_gra)):
        f.write('{} {} {} {}\n'.format(yp[i],xp[i],zp[i],np.array(field_gra[i]).ravel()[0]))
        
with open(magoutfile, 'w') as f:
    f.write('! model 2 magtotal-field magnetic anomaly (nT)\n')
    f.write('{}\n'.format(len(field_mag)))
    for i in range(len(field_mag)):
        f.write('{} {} {} {}\n'.format(yp[i],xp[i],zp[i],np.array(field_mag[i]).ravel()[0]))
    
with open(graoutfile1, 'w') as f:
    f.write('! model 2 gravity anomlay (mGal) with 5% noise\n')
    f.write('{}\n'.format(len(field_gra1)))
    for i in range(len(field_gra1)):
        f.write('{} {} {} {}\n'.format(yp[i],xp[i],zp[i],np.array(field_gra1[i]).ravel()[0]))
    
with open(magoutfile1, 'w') as f:
    f.write('! model 2 magtotal-field magnetic anomaly (nT) with 5% noise\n')
    f.write('{}\n'.format(len(field_mag1)))
    for i in range(len(field_mag1)):
        f.write('{} {} {} {}\n'.format(yp[i],xp[i],zp[i],np.array(field_mag1[i]).ravel()[0]))
    
    
#画图
plt.figure(figsize=(16, 16))
plt.subplot(2, 2, 1)
plt.axis('scaled')
plt.title('model 2 gravity anomlay (mGal)')
levels = giplt.contourf(yp , xp , field_gra, nshape, 15)
cb = plt.colorbar(orientation='horizontal')
giplt.contour(yp, xp, field_gra, nshape,
                levels, clabel=False, linewidth=0.1)
plt.subplot(2, 2, 2)
plt.axis('scaled')
plt.title('model 2 magtotal-field magnetic anomaly (nT)')
levels = giplt.contourf(yp , xp , field_mag, nshape, 15)
cb = plt.colorbar(orientation='horizontal')
giplt.contour(yp, xp, field_mag, nshape,
                levels, clabel=False, linewidth=0.1)

plt.subplot(2, 2, 3)
plt.axis('scaled')
plt.title('model 2 gravity anomlay (mGal) with 5% noise')
levels = giplt.contourf(yp , xp , field_gra1, nshape, 15)
cb = plt.colorbar(orientation='horizontal')
giplt.contour(yp, xp, field_gra1, nshape,
                levels, clabel=False, linewidth=0.1)
plt.subplot(2, 2, 4)
plt.axis('scaled')
plt.title('model 2 magtotal-field magnetic anomaly (nT) with 5% noise')
levels = giplt.contourf(yp , xp , field_mag1, nshape, 15)
cb = plt.colorbar(orientation='horizontal')
giplt.contour(yp, xp, field_mag1, nshape,
                levels, clabel=False, linewidth=0.1)
plt.tight_layout()
plt.show()