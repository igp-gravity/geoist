# -*- coding: utf-8 -*-
"""
Created on Sun Mar 10 20:49:28 2019

@author: chens
"""
import numpy as np
from scipy import sparse

def haarMatrix(n, normalized=False):
    # Allow only size n of power 2
    n = 2**np.ceil(np.log2(n))
    if n > 2:
        h = haarMatrix(n / 2)
    else:
        if normalized:
            return np.dot(np.sqrt(2)/2, np.array([[1, 1], [1, -1]]))
        else:
            return np.array([[1, 1], [1, -1]])       
    # calculate upper haar part
    h_n = np.kron(h, [1, 1])
    # calculate lower haar part 
    if normalized:
        h_i = np.sqrt(n/2)*np.kron(np.eye(len(h)), [1, -1])
        h = np.dot(1/np.sqrt(n/2),np.vstack((h_n, h_i)))
    else:
        h_i = np.kron(np.eye(len(h)), [1, -1])
        h = np.vstack((h_n, h_i))
    return h

def walshMatrix(n, normalized=False):
    # Allow only size n of power 2
    n = 2**np.ceil(np.log2(n))
    if n > 2:
        h = walshMatrix(n / 2)
    else:
        if normalized:
            return np.dot(1/np.sqrt(2), np.array([[1, 1], [1, -1]]))
        else:
            return np.array([[1, 1], [1, -1]])
    # calculate upper haar part
    h_n = np.kron([[1, 1], [1, -1]], h)
    # calculate lower haar part 
    if normalized:
        h = np.sqrt(2)**-np.log2(n)*h_n
    else:
        h = h_n
    return h

if __name__ == '__main__':

    import time
    #import pywt
    nn = 8
    time_start=time.time()
    a1 = np.mat(walshMatrix(nn, normalized=False))
    #np.linalg.svd(a1)
    #np.linalg.det(a1)
    print(a1)
    time_end=time.time()
    print('totally cost',time_end-time_start)
#    aa1 = np.where(np.abs(a1*a1.T) > 1.0e-10, a1*a1.T, 0)
#    print(aa1)    
#    ax = np.mat(walshMatrix(nn, normalized=True))
#    aax = np.where(np.abs(ax*ax) > 1.0e-10, ax*ax, 0)
#    #print(aax)
#    a=np.array([(10,9,6,5),(9,10,9,6),(6,9,10,9),(5,6,9,10)])
#    h1=np.sqrt(2)*np.array([(0.5,0.5,0.0,0.0),
#                (0.0,0.0,0.5,0.5),
#                (0.5,-0.5,0.0,0.0),
#                (0.0,0.0,0.5,-0.5)])
    h1=np.array([(1,1,0,0),
                (0,0,1,1),
                (1,-1,0,0),
                (0,0,1,-1)])
#    a=np.mat(a)
#    sA = sparse.csr_matrix(a)
#    h1=np.mat(h1)
#    c0=np.dot(a,h1.T)
#    d0=np.dot(h1,c0)
#    c1=np.dot(d0,h1.T)
#    d1=np.dot(h1,c1)
#    d1=h1*h1*a*h1.T*h1.T
##    print(h1*h1*a*h1.T*h1.T)   
##    print(len(np.nonzero(d1)[0]))
#    
##    ggz = np.loadtxt('d:\\ggz888.dat')
##    gz = np.zeros(np.shape(ggz))
##    j = 0
##    for i in range(int(nn/2)):
##        gz[i] = ggz[j]
##        gz[i+int(nn/2)] = ggz[j+1]
##        j = j + 2  
##    gz = np.mat(ggz)
##    d2=ax*gz.T*gz*ax.T
##    
##    #d2=ax*ax*ax*ax*ax*ax*ax*ax*ax*gz.T*gz*ax.T*ax.T*ax.T*ax.T*ax.T*ax.T*ax.T*ax.T*ax.T
##    d3=np.where(np.abs(d2) > 1.0e-10, d2, 0)
##    print(len(np.nonzero(d3)[0]))
##    a2=np.hstack((a,0.5*a))
##    a4=a2.T*a2 #np.vstack((a2,-0.3*a2))
##    a4=np.mat(a4)
#    h2=np.sqrt(2)*np.array([(0.5,0.5,0.0,0.0,0.0,0.0,0.0,0.0),
#                 (0.0,0.0,0.5,0.5,0.0,0.0,0.0,0.0),
#                 (0.0,0.0,0.0,0.0,0.5,0.5,0.0,0.0),
#                 (0.0,0.0,0.0,0.0,0.0,0.0,0.5,0.5),
#                 (0.5,-0.5,0.0,0.0,0.0,0.0,0.0,0.0),
#                 (0.0,0.0,0.5,-0.5,0.0,0.0,0.0,0.0),
#                 (0.0,0.0,0.0,0.0,0.5,-0.5,0.0,0.0),
#                 (0.0,0.0,0.0,0.0,0.0,0.0,0.5,-0.5)]) 
#    h2=np.mat(h2)
##    d2=h2*h2*h2*a4*h2.T*h2.T*h2.T
##    d3=np.where(np.abs(d2) > 1.0e-10, d2, 0)
##    print(d3)
##    print(len(np.nonzero(d3)[0]))
##    b=np.array([(88,88,89,90,92,94,96,97),
##                (90,90,91,92,93,95,97,97),
##                (92,92,93,94,95,96,97,97),
##                (93,93,94,95,96,96,96,96),
##                (92,93,95,96,96,96,96,95),
##                (92,94,96,98,99,99,98,97),
##                (94,96,99,101,103,103,102,101),
##                (95,97,101,104,106,106,105,105)])
##    coeffs1 = pywt.wavedec2(b, 'db1','periodization')
##    pywt.waverec2(coeffs1, 'db1')
##    xx,ss = pywt.coeffs_to_array(coeffs1)
##    print(xx)