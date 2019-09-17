# -*- coding: utf-8 -*-
"""
Created on Sun Mar 10 20:49:28 2019

@author: chens
"""
import numpy as np

def walshMatrix(n, normalized=False):
    # Allow only size n of power 2
    n = 2**np.ceil(np.log2(n))
    n = int(n)
    h = np.zeros((n,n))
    j = 0
    for i in range(n):
        if i <= n/2 - 1:
            h[i,j] = 0.5
            h[i,j+1] = 0.5
            j = j + 2
        else:
            if i == n/2:
                j = 0
            h[i,j] = 0.5
            h[i,j+1] = -0.5 
            j = j + 2
    if normalized:        
        h = np.mat(np.sqrt(2)*h)
    else:
        h = np.mat(h)
    h0 = np.mat(np.eye(n))
    for k in range(int(np.log2(n))):
        h0 = h0*h
    return h0, h

if __name__ == '__main__':
    from scipy import sparse
    nn = 512
    ax, ax0 = walshMatrix(nn, normalized=True)
    sax0 = sparse.csr_matrix(ax0)
    #ax, ax0 = np.mat(walshMatrix(nn, normalized=True))
    
#    wm = []
#    for k in range(1,9):
#        wm = np.hstack((wm,k*np.ones(64)))        
#    wm = np.mat(np.diag(wm))
#    mref = wm*np.eye(nn, dtype=int)
#    ggz0 = np.mat(np.loadtxt('d:\\ggz888.dat'))
#    ggz = np.vstack((ggz0,mref))
#    d2 = ax*ggz.T*ggz*ax.T
#    gz0 = np.mat(np.loadtxt('d:\\gz888.dat')).T
#    obsref = np.mat(np.zeros(nn))
#    gz = np.vstack((gz0,obsref.T))
#    gzw = ax*ggz.T*gz
#    d3 = np.where(np.abs(d2) > 1.0e-10, d2, 0)
#    print(len(np.nonzero(d3)[0])) #非零元素
#    if np.allclose(d3, d3.T, atol=1e-8):
#        print('压缩后的矩阵对称')
#    g4 = np.where(np.abs(gzw) > 1.0e-10, gzw, 0)
#    print(len(np.nonzero(g4)[0]))
#    m1 = np.linalg.solve(d3, g4)
#    m2 = ax*m1
#    np.savetxt("d:\\walshinv\\m2.txt", m2)
#    if np.allclose(ggz0*m2, gz0, atol=1e-2):
#        print('Result is right')
#        x1 = ggz0*m2-gz0
#        print(np.sqrt(x1.T*x1))
#    else:
#        x1 = ggz0*m2-gz0
#        print(np.sqrt(x1.T*x1))
#    m3 = np.linalg.lstsq(ggz, gz, rcond=None)
#    np.savetxt("d:\\walshinv\\m3.txt", m3[0])
#    print(np.allclose(ggz*m3[0], gz)) 