# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 12:50:51 2019
   Small Angle Approximation
@author: chens
"""
import numpy as np

alfa = 0.5
beta = 0.2
xita = -0.3

alfa = alfa*np.pi/180. #Z
beta = beta*np.pi/180. #Y
xita = xita*np.pi/180. #X
# small angle approximation
a1 = 1.
a2 = -alfa
a3 = beta
a4 = alfa
a5 = 1.
a6 = -xita
a7 = -beta
a8 = xita
a9 = 1.
rmat = np.mat([[a1, a2, a3],[a4, a5, a6],[a7, a8, a9]])

# rotation matrices
t1 = np.cos(beta)*np.cos(alfa) + np.sin(beta)*np.sin(alfa)*np.sin(xita)
t2 = -np.cos(alfa)*np.sin(alfa) + np.sin(beta)*np.cos(alfa)*np.sin(xita)
t3 = np.sin(beta)*np.cos(xita)
t4 = np.sin(alfa)*np.cos(xita)
t5 = np.cos(alfa)*np.cos(xita)
t6 = -np.sin(xita)
t7 = -np.sin(beta)*np.cos(alfa) + np.cos(beta)*np.sin(alfa)*np.sin(xita)
t8 = np.sin(beta)*np.sin(alfa) + np.cos(beta)*np.cos(alfa)*np.sin(xita)
t9 = np.cos(beta)*np.cos(xita)
tmat = np.mat([[t1, t2, t3],[t4, t5, t6],[t7, t8, t9]])


# there are must two points in the transform.
p1 = np.array((027747.74, 001439.14, 036583.48))
p11 = np.array((027747.33, 001438.91, 036584.64))
# the coordinates in the new references
p2 = np.array((027075.90, 001436.96, 037131.64)) #rmat*np.transpose(np.mat(p1))
p21 = np.array((027075.91, 001437.01, 037131.94)) #rmat*np.transpose(np.mat(p11))
#print(p2)
mcore = np.zeros(54).reshape(6,9)
mcore[0, 0:3] = p1
mcore[1, 3:6] = p1
mcore[2, 6:9] = p1
mcore[3, 0:3] = p11
mcore[4, 3:6] = p11
mcore[5, 6:9] = p11
mext = np.zeros(54).reshape(6,9)
mext[0, 0] = 1
mext[1, 4] = 1
mext[2, 8] = 1
mext[3, 1] = 1
mext[3, 3] = 1
mext[4, 2] = 1
mext[4, 6] = 1
mext[5, 5] = 1
mext[5, 7] = 1
mmker = np.vstack((mcore, mext))
#print(mmker) 
pobs = np.zeros(12)
pobs[0:3] = p2.T
pobs[3:6] = p21.T
pobs[6:9] = 1.0
#print(p2)
print(np.linalg.matrix_rank(mmker.T@mmker))
x0 = np.linalg.solve(mmker.T@mmker, mmker.T@pobs)
print(x0)

#动态情况, Δt时间间隔除了转动还有位置变化，即两个时间的M矢量距离不相等
bb = np.diag(np.ones(3))  #增加磁异常梯度可以测量或者读出来
bbb = np.vstack((bb, bb))
xx1 = np.hstack((mcore, bbb))
xx2 = np.hstack((mext, np.zeros(18).reshape(6,3)))
xx3 = np.vstack((xx1, xx2))
np.linalg.matrix_rank(xx3) #不满阵，需要加入加速度传感器约束

xx4 = np.zeros(36).reshape(3,12)  #增加加速度传感器积分后的信息
xx4[0,2] = 5
xx4[1,3] = 5
xx4[2,7] = 5   #假设5为加速度传感器积分得到的位移矢量
xx4[0,9] = -1
xx4[1,10] = -1
xx4[2,11] = -1

xx4 = np.zeros(36).reshape(3,12)
xx5 = np.vstack((xx3, xx4))  #得到动态条件下的小角度坐标平差方程
np.linalg.matrix_rank(xx5)  
