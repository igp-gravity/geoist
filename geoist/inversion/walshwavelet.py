# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 11:06:23 2019

@author: chens
"""

import numpy as np
from scipy import sparse

def grav_matrix(number, normalized=False):
    '''
    generate walsh matrix
    
    '''
    # Allow only size n of power 2
    number = 2**np.ceil(np.log2(number))
    number = int(number)
    h_w = sparse.lil_matrix((number, number), dtype=np.int8)
    j = 0
    for i in range(number):
        if i <= number/2 - 1:
            h_w[i, j] = 1
            h_w[i, j+1] = 1
            j = j + 2
        else:
            if i == number/2:
                j = 0
            h_w[i, j] = 1
            h_w[i, j+1] = -1
            j = j + 2
#    if normalized:
#        h = np.mat(np.sqrt(2)*h)
#    else:
#        h = np.mat(h)
#    h0 = np.mat(np.eye(n))
#    for k in range(int(np.log2(n))):
#        h0 = h0*h
    #sh = sparse.csr_matrix(h)
    return h_w

if __name__ == '__main__':
    import sys
    NN = 64*64*32
    #h = sparse.lil_matrix((nn, nn))
    AX = grav_matrix(NN, normalized=False)
    print(NN)
    print(sys.getsizeof(AX))
