# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 10:40:28 2019

@author: chens
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv(r'D:\\AbsGra\\20190619.dat',sep=' ',encoding='utf-8')     
mag_obs_3 = df.iloc[:,[4,5,3]]
atemp = np.array(mag_obs_3[0:79]).astype(np.float64)
magabs = np.sqrt(atemp[:,0]**2+atemp[:,1]**2+atemp[:,2]**2)
#plt.plot(magabs)
#plt.show()

# pitch，roll，yaw/azimuth/heading / angle in radian  y x z/ fai omega kai
# 2.08 2.05 2.28
pitch = 1.17*np.pi/180. # Y
roll = 0.82*np.pi/180. # X 
yaw = 1.07*np.pi/180. # Z

# Magnetic X/North Y/East Z/Vertical 
Mag_x, Mag_y, Mag_z = ( 027747.74, 001439.14, 036583.48 )
Mag_x1, Mag_y1, Mag_z1 = ( 027075.90, 001436.96, 037131.64 )

print(Mag_x, Mag_y, Mag_z)
print(Mag_x1, Mag_y1, Mag_z1)
# small angle approximation
a1 = 1.
a2 = yaw  
a3 = -pitch
a4 = -yaw
a5 = 1.
a6 = roll
a7 = pitch
a8 = -roll
a9 = 1.
rmat = np.mat([[a1, a2, a3],[a4, a5, a6],[a7, a8, a9]])
# rotation matrices
t1 = np.cos(pitch)*np.cos(yaw)
t2 = np.cos(roll)*np.sin(yaw) + np.sin(roll)*np.sin(pitch)*np.cos(yaw)
t3 = np.sin(roll)*np.sin(yaw) - np.cos(roll)*np.sin(pitch)*np.cos(yaw)
t4 = -np.cos(pitch)*np.sin(yaw)
t5 = np.cos(roll)*np.cos(yaw) - np.sin(roll)*np.sin(pitch)*np.sin(yaw)
t6 = np.sin(roll)*np.cos(yaw) + np.cos(roll)*np.sin(pitch)*np.sin(yaw)
t7 = np.sin(pitch)
t8 = -np.sin(roll)*np.cos(pitch)
t9 = np.cos(roll)*np.cos(pitch)

tmat = np.mat([[t1, t2, t3],[t4, t5, t6],[t7, t8, t9]])
# there are must two points in the transform.
p1 = np.array(( Mag_x, Mag_y, Mag_z))
# the coordinates in the new references
p2 = tmat*np.transpose(np.mat(p1))
print(np.array(p2.reshape(3)))
print(np.sqrt(Mag_x**2 + Mag_y**2 + Mag_z**2))
print(np.sqrt(p2[0]**2 +  p2[1]**2 +  p2[2]**2))