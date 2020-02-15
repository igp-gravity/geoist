import numpy as np
from numba import jit

@jit(nopython=True,nogil=True,parallel=True)
def natural_walsh_matrix_old(n,normalized=False):
    '''natural ordering walsh matrix

    Args:
        n (int): degree of walsh matrix.
        normalized (bool): whether or not normalize the generated matrix.

    Returns:
        walsh matrix in natural ordering.
    '''
    n = 2**np.ceil(np.log2(n))
    n = int(n)
    h = np.zeros((n,n))
    j = 0
    for i in range(n):
        if i <= n/2 -1:
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
        h = np.sqrt(2)*h
    h0 = np.eye(n)
    for k in range(int(np.log2(n))):
        h0 = np.dot(h0,h)
    return h0

#@jit(nopython=True,nogil=True,parallel=True)
def natural_walsh_matrix(n,normalized=False):
    '''natural ordering walsh matrix

    Args:
        n (int): degree of walsh matrix.
        normalized (bool): whether or not normalize the generated matrix.

    Returns:
        walsh matrix in natural ordering.
    '''
    logn = int(np.ceil(np.log2(n)))
    n = 2**logn
    h = np.array([[1.,1.],[1.,-1.]])
    for i in range(logn-1):
        h = np.vstack((np.hstack((h,h)),np.hstack((h,-h))))
    if normalized:
        h = h/(np.sqrt(2)**logn)
    return h

def walsh_order(n):
    ''' generate 'natural','dyadic','sequence' ordering of walsh matrix.

    Args:
        n (int): degree of walsh matrix.
    '''
    n = 2**np.ceil(np.log2(n))
    n = int(n)
    n_bits = len(np.binary_repr(n))-1
    print(n_bits)
    sequence_order = np.arange(n)
    tmp = np.right_shift(sequence_order,1)
    dyadic_order = np.bitwise_xor(sequence_order,tmp)
    natural_order = [int('{:0{width}b}'.format(i,width=n_bits)[::-1],2) for i in dyadic_order]
    return sequence_order,dyadic_order,natural_order

def walsh_matrix_reorder(m,origin_order='natural',target_order='sequence'):
    ''' reorder walsh matrix

    Args:
        m (ndarray): walsh matrix
        from (string): orginal ordering. take values from 'natural','dyadic',
            'sequence'
        to (string): target ordering. take values from 'natural','dyadic',
            'sequence'

    Returns:
        reordered walsh matrix.
    '''
    sequence_order,dyadic_order,natural_order = walsh_order(m.shape[0])
    order_dict = {'natural':natural_order,'dyadic':dyadic_order,'sequence':sequence_order}
    tmp_order = sequence_order.copy()
    tmp_order[order_dict[target_order]] = order_dict[origin_order]
    m[:] = m[tmp_order]

def walsh_reorder2(m,nxyz):
    ''' reorder walsh matrix from sequence order to sequence2 order. In
    sequence2 ordering, kernel.T*kernel will be a block diagonal matrix.

    Args:
        m (ndarray): walsh matrix in sequence ordering.
        nxyz (tuple of ints): (nx,ny,nz) specifies source points configuration.

    Returns:
        walsh matrix in sequence2.

    '''
    nx,ny,nz = nxyz
    big_block_n = nz*ny
    small_block_n = nz
    # first reorder big blocks of m.
    big_pre = np.arange(0,nx,2)*big_block_n # non-zero block start index
    big_post = np.arange(1,nx,2)*big_block_n # zero-block start index
    big_order_1 = np.hstack(tuple(np.arange(big_pre[i],big_post[i]) for i in range(len(big_pre))))
    big_pre += 2*big_block_n
    big_order_2 = np.hstack(tuple(np.arange(big_post[i],big_pre[i]) for i in range(len(big_pre))))
    big_order = np.hstack((big_order_1,big_order_2))
    # then reorder small blocks inside each big block.
    small_pre = np.arange(0,ny*nx//2,2)*small_block_n
    small_post = np.arange(1,ny*nx//2,2)*small_block_n
    small_order_1 = np.hstack(tuple(np.arange(small_pre[i],small_post[i]) for i in range(len(small_pre))))
    small_pre += 2*small_block_n
    small_order_2 = np.hstack(tuple(np.arange(small_post[i],small_pre[i]) for i in range(len(small_pre))))
    small_order = np.hstack((small_order_1,small_order_2,small_order_1+nx*ny*nz/2,small_order_2+nx*ny*nz/2)).astype(int)
    m[:] = m[big_order][small_order]

def walsh_matrix(n,normalized=False,ordering='sequence',nxyz=(None,None,None)):
    ''' generate walsh matrix with given ordering.

    Args:
        n (int): degree of walsh matrix.
        normalized (bool): whether or not nomalize walsh matrix.
        ordering (string): ordering of walsh matrix. take values from 'natural',
            'dyadic','sequence'
        nxyz (tuple of ints): (nx,ny,nz) specifies source points configuration.

    Returns:
        walsh matrix with given ordering.

    '''
    m = natural_walsh_matrix(n,normalized)
    if ordering != 'sequence2':
        walsh_matrix_reorder(m,target_order=ordering)
    else:
        walsh_matrix_reorder(m,target_order='sequence')
        walsh_reorder2(m,nxyz)
    return m

@jit(nopython=True,nogil=True,parallel=True)
def linear_transform(M,IM,A=None,x=None,b=None):
    ''' Transform equation Ax=b to A'x'=b'
    where A'=M*A*M^{-1}
          x'=M*x
          b'=M*b.

    Args:
        M (ndarray): trasform matrix.
        IM (ndarray): inverse of trasform matrix.
        A (ndarray): coefficient matrix.
        x (ndarray): unkowns.
        b (ndarray): righthand side.

    Returns:
        transformed A',x',b'.
    '''
    if not A is None:
        A = M @ A @ IM
    if not x is None:
        x = M @ x
    if not b is None:
        b = M @ b
    return A,x,b

@jit(nopython=True,nogil=True,parallel=True)
def linear_recover(M,IM,A=None,x=None,b=None):
    ''' Recover equation Ax=b from A'x'=b'
    where A'=M*A*M^{-1}
          x'=M*x
          b'=M*b.

    Args:
        M (ndarray): trasform matrix.
        IM (ndarray): inverse of trasform matrix.
        A (ndarray): coefficient matrix A'.
        x (ndarray): unkowns x'.
        b (ndarray): righthand side b'.

    Returns:
        recovered A,x,b.
    '''
    if not A is None:
        A = IM @ A @ M
    if not x is None:
        x = IM @ x
    if not b is None:
        b = IM @ b
    return A,x,b

@jit(nopython=True,nogil=True,parallel=True)
def walsh_transform(A=None,x=None,b=None,normalized=True,
                    ordering='sequence2',nxyz=(None,None,None)):
    ''' Transform linear equation Ax=b to A'x'=b'
    where A'=M*A*M^{-1}
          x'=M*x
          b'=M*b.
          M is the walsh matrix.

    Args:
        n (int): degree of walsh matrix.
        M (ndarray): trasform matrix.
        A (ndarray): coefficient matrix.
        x (ndarray): unkowns.
        b (ndarray): righthand side.
        ordering (string): ordering of walsh matrix.
          could be one of 'natural','dyadic','sequence','sequence2'.
        nxyz (tuple of ints): (nx,ny,nz) specifies source points configuration.

    Returns:
        transformed A',x',b'.
    '''
    if M is None:
        M = walsh_matrix(n,normalized=normalized,ordering=ordering,nxyz=nxyz)
    return linear_transform(M.T,M,A,x,b)

@jit(nopython=True,nogil=True,parallel=True)
def walsh_recover(n,M=None,A=None,x=None,b=None,normalized=True,
                    ordering='sequence2',nxyz=(None,None,None)):
    ''' Recover linear equation Ax=b from A'x'=b'
    where A'=M*A*M^{-1}
          x'=M*x
          b'=M*b.
          M is the walsh matrix.

    Args:
        n (int): degree of walsh matrix.
        M (ndarray): trasform matrix.
        A (ndarray): coefficient matrix.
        x (ndarray): unkowns.
        b (ndarray): righthand side.
        ordering (string): ordering of walsh matrix.
          could be one of 'natural','dyadic','sequence','sequence2'.
        nxyz (tuple of ints): (nx,ny,nz) specifies source points configuration.

    Returns:
        recovered A,x,b.
    '''
    if M is None:
        M = walsh_matrix(n,normalized=normalized,ordering=ordering,nxyz=nxyz)
    return linear_recover(M.T,M,A,x,b)

@jit(nopython=True,nogil=True,parallel=True)
def fast_walsh_transform(a,normalized=True):
    '''
    inplace FFT-like walsh transform. Let M be the walsh matrix, this routine
    will return M*a. This routine is designed to transform a vector. Although it
    still works on matrix, it misses the multiplier M^{-1} on the right side. If
    you want transform matrix A, first run fast_walsh_transform(A), then
    fast_walsh_transform(A.T).
    This routine use natural ordering walsh matrix.

    Args:
        a (ndarray): ndarray to be transformed

    Returns:
        M*a. where M is the walsh matrix
    '''
    h = 1
    if normalized:
        coef = np.sqrt(2)
    else:
        coef = 1.
    while h < len(a):
        for i in range(0,len(a),h*2):
            for j in range(i,i+h):
                x = a[j]
                y = a[j+h]
                a[j] = (x + y) / coef
                a[j+h] = (x - y) / coef
        h *= 2
