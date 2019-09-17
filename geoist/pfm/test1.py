# -*- coding: utf-8 -*-
"""
Created on Wed May  1 18:44:52 2019

@author: chens
"""
from scipy.linalg import toeplitz
import numpy as np
t1 = toeplitz([1,2,3,4])
t2 = toeplitz([6,7,3,9])
#scipy.sparse.linalg.eigsh #Find k eigenvalues of the real symmetric square matrix 
h1=np.array([(1,1,0,0),
            (0,0,1,1),
            (1,-1,0,0),
            (0,0,1,-1)])
h2=np.array([(1,1,0,0,0,0,0,0),
            (0,0,1,1,0,0,0,0),
            (0,0,0,0,1,1,0,0),
            (0,0,0,0,0,0,1,1),
            (1,-1,0,0,0,0,0,0),
            (0,0,1,-1,0,0,0,0),
            (0,0,0,0,1,-1,0,0),
            (0,0,0,0,0,0,1,-1)])
h1 = np.mat(h1)
h2 = np.mat(h2)
w4 = h1*h1
w8 = h2*h2*h2
tt = np.hstack((t1,t2))
tt = np.mat(tt)
hh = np.vstack((h1,h1))
b1=np.mat([(1,-1,0,0),
            (0,1,-1,0),
            (0,0,1,-1)])
b2=np.mat([(1,-1,0,0),
            (0,1,-1,0),
            (0,0,1,-1)])
bb = np.hstack((b1,b2))


X = np.random.rand(4**2).reshape(4, 4)
X = np.triu(X)
X += X.T - np.diag(X.diagonal())