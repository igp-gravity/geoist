import numpy as np
from scipy import linalg as splin
import numba as nb
from scipy import fftpack as spfft
#import cupy
import time
use_gpu = 0

# if use_gpu > 0:
#     import cupy
# else:
#     print('NO GPU BE USED',use_gpu)

def check_gpu():
    print(use_gpu)
    
def block_circ(a):
    '''generate full representation of 2-level circulant matrix

    Args:
        a (ndarray): 1-st column of the circulant matrix in proper shape.
    Returns:
        Full filled circulant matrix
    '''
    if use_gpu > 0:
        import cupy
        xp = cupy.get_array_module(a)
        if xp is cupy:
            a = xp.asnumpy(a)
    else:
        xp = np
        a = np.asnumpy(a)
            
    if a.ndim == 1:
        return splin.circulant(a)
    n_total = np.prod(a.shape)
    a_shape = a.shape
    A = np.zeros(np.hstack([np.array(a.shape),np.array(a.shape)]))
    A_shape = A.shape
    x_slice = [0]*a.ndim
    x_slice[-1] = slice(None)
    x_target_slice = [0]*a.ndim
    x_target_slice[-1] = slice(None)
    y_slice = [slice(None)]*2
    a = a.reshape(-1,a_shape[-1])
    A = A.reshape(-1,a_shape[-1],*a_shape)
    for i,sub_column in enumerate(a):
        print(sub_column)
        y_slice[0] = i
        A[tuple(y_slice+x_slice)] = splin.circulant(sub_column)
    for ilevel in range(len(a_shape)-1):
        A = A.reshape(-1,*a_shape[len(a_shape)-ilevel-2:],*a_shape)
        y_slice = [slice(None)]*(ilevel+3)
        y_target_slice = [slice(None)]*(ilevel+3)
        for k in range(A.shape[0]):
            y_slice[0] = k
            y_target_slice[0] = k
            for i in range(a_shape[len(a_shape)-ilevel-2]):
                y_slice[1] = i
                for j in range(1,a_shape[len(a_shape)-ilevel-2]):
                    x_slice[len(a_shape)-ilevel-2] = j
                    y_target_slice[1] = np.mod(i-j,a_shape[len(a_shape)-ilevel-2])
                    A[tuple(y_slice+x_slice)] = A[tuple(y_target_slice+x_target_slice)]
        x_slice[len(a_shape)-ilevel-2] = slice(None)
        x_target_slice[len(a_shape)-ilevel-2] =slice(None)
    A = A.reshape(A_shape)
    return A

def block_toep2_sym(a):
    '''generate full representation of 2-level symmetric toeplitz matrix

    Args:
        a (ndarray): 1-st column of the symmetrec toeplitz matrix in proper shape.
    Returns:
        Full filled toeplitz matrix.
    '''
    if use_gpu > 0:
        import cupy
        xp = cupy.get_array_module(a)
        if xp is cupy:
            a = xp.asnumpy(a)
    else:
        xp = np
        a = np.asnumpy(a)
        
    A1 = []
    n0,n1 = a.shape
    for i in range(n1):
        A1.append(splin.toeplitz(a[:,i]))
    A = np.empty((n1,n0,n1,n0))
    for i in range(n1):
        for j in range(n1):
            A[i,:,j,:] = A1[np.int(np.abs(i-j))]
    A.shape = (n0*n1,n0*n1)
    A = xp.asarray(A)
    
    return(A)

def circ_eigs(circ,dtype=np.complex64):
    ''' calculate eigenvalues of multilevel circulant matrix A

    Args:
        circ (ndarray): representation of the multilevel circulant matrix A, i.e.
             the first column of A in proper shape.
    Returns:
        eigenvalues of A.
    '''
    if use_gpu > 0:
        import cupy
        xp = cupy.get_array_module(circ)
    else:
        xp = np
    return xp.fft.fft2(circ,norm='ortho').astype(dtype)*xp.sqrt(np.prod(circ.shape))

def circ_mul_v(circ,v,eigs=None):
    ''' multiply a circulant matrix A by multi vector v.

    Args:
        circ (ndarray): representation of the multilevel circulant matrix A, i.e.
             the first column of A in proper shape.
        v (ndarray): vector to be multiplied. Should be reshaped to the same shape
             as circ. Should be the same reshape order (column first/row first) as circ.
    Returns:
        result of multiplication.
    '''
    if use_gpu > 0:
        import cupy
        xp = cupy.get_array_module(circ)
    else:
        xp = np
    
    if eigs is None:
        eigs = circ_eigs(circ)
    tmp = xp.real(xp.fft.ifft2(xp.fft.fft2(v,norm='ortho')*eigs,norm='ortho'))
    if xp is cupy:
        return tmp.astype(xp.float32)
    else:
        return tmp

def embed_toep2circ(toep,v=None):
    '''embed toeplitz matrix to circulant matrix.

    Args:
        toep (ndarray): representation of multilevel toeplitz matrix, i.e.
            the first column of A in proper shape.
        v (ndarray): embed a vector together with toep.
    Returns:
        representation of embedded multilevel circulant matrix and embedded vector.
    '''
    if use_gpu > 0:
        import cupy
        xp = cupy.get_array_module(toep)
    else:
        xp = np
        
#    circ = xp.zeros((2*xp.array(toep.shape)).astype(xp.int))
    circ = xp.zeros((2*toep.shape[0],2*toep.shape[1]))
    s = []
    for idim in range(toep.ndim):
        s.append(slice(0,toep.shape[idim]))
    circ[tuple(s)] = toep
    if not v is None:
        if v.ndim == toep.ndim:
            resv = xp.zeros_like(circ)
            resv[tuple(s)] = v
        else:
            resv = xp.zeros((v.shape[0],*circ.shape))
            resv[tuple(slice(None),*s)] = v
    for idim in range(toep.ndim):
        s[idim] = slice(toep.shape[idim]+1,None)
        s2 = s.copy()
        s2[idim] = slice(toep.shape[idim]-1,0,-1)
        circ[tuple(s)] = circ[tuple(s2)]
        s[idim] = slice(None)
    if v is None:
        return circ
    else:
        return circ,resv

def toeplitz_mul_v(toep,v):
    '''multiply a toeplitz matrix A by vector v.

    Args:
        toep (ndarray): representation fo multilevel toeplitz matrix, i.e.
            the first column of A in proper shape.
        v (ndarray): vector to be multiplied. Should be reshaped to the same shape
             as toep. Should be the same reshape order (column first/row first) as circ.
    Returns:
        result of multiplication.
    '''
    
    #xp = cupy.get_array_module(toep)
    circ,tmpv = embed_toep2circ(toep,v)
    res = circ_mul_v(circ,tmpv)
    slices = [slice(None)]*res.ndim
    slices[-1] = slice(0,v.shape[-1])
    slices[-2] = slice(0,v.shape[-2])
    return(res[tuple(slices)])

class GToepOperator:
    '''This class construct an linear operator of multilevel symmetric Toeplitz
       matrix.'''

    def __init__(self,toep,*args,**kwargs):
        '''Input toep is the first row of the underlying Toeplitz matrix'''
        if use_gpu > 0:
            import cupy
            self.xp = cupy.get_array_module(toep)
        else:
            self.xp = np
            
        self.nz,self.ny,self.nx = toep.shape
        self.eigs = self.xp.zeros((self.nz,2*self.ny,2*self.nx),dtype=np.complex64)
        for i in range(self.nz):
            tmptoep = embed_toep2circ(toep[i])
            self.eigs[i] = circ_eigs(tmptoep,dtype=np.complex64)
        self.dtype = toep.dtype
        self.get_m_eigs(toep)

    def get_m_eigs(self,toep):
        '''M is designed to be an matrix approximate the underlying Toeplitz
           matrix. So that M could be used as preconditioner.'''
        M = self.xp.zeros((self.ny,self.nx))
        for i in range(self.nz):
            M += toep[i]
        for idim in range(M.ndim):
            coef = (1.*M.shape[idim]-self.xp.arange(1,M.shape[idim]))/M.shape[idim]
            coef_shape = np.ones(M.ndim,dtype=np.int)
            coef_shape[idim] = M.shape[idim] - 1
            coef = coef.reshape(tuple(coef_shape))
            s1 = [slice(None)]*M.ndim
            s2 = [slice(None)]*M.ndim
            s1[idim] = slice(1,None)
            s2[idim] = slice(-1,None,-1)
            s1 = tuple(s1)
            s2 = tuple(s2)
            M[s1] = M[s1]*coef + (M[s1]*coef)[s2]
        self.m_eigs = circ_eigs(M)
#        self.m_eigs = self.xp.where(self.xp.abs(self.m_eigs)>1./self.xp.sqrt(self.ny*self.nx),
#                                    self.m_eigs,
#                                    (1.+0.j)/self.xp.sqrt(self.ny*self.nx))

    def matvec_prec(self,v):
        ''' Apply M^{-1} to a vector v.'''
        if v.ndim == 1:
            return circ_mul_v(self.m_eigs,v.reshape(self.ny,self.nx),1./self.m_eigs)
        else:
            return circ_mul_v(self.m_eigs,v.reshape(-1,self.ny,self.nx),1./self.m_eigs)

    def matvec_prec_sym(self,v):
        ''' Apply M^{-T}M^{-1} to a vector v.'''
        if use_gpu > 0:
            import cupy
        tmp = self.xp.fft.fft2(v.reshape(-1,self.ny,self.nx),norm='ortho')
        tmp = tmp/(self.m_eigs*self.xp.conj(self.m_eigs))
        tmp = self.xp.real(self.xp.fft.ifft2(tmp,norm='ortho'))
        if v.ndim == 1:
            tmp = tmp.ravel()
        else:
            tmp = tmp.reshape(v.shape[0],-1)
        if self.xp is cupy:
            return tmp.astype(self.xp.float32)
        else:
            return tmp

    def matvec(self,v):
        ''' Apply the underlying Toeplitz matrix T to a multivector v.'''
        v_orig = v.reshape(-1,self.nz,self.ny,self.nx)
        res = self.xp.zeros((v_orig.shape[0],2*self.ny,2*self.nx),dtype=np.complex64)
        expand_v = np.zeros((v_orig.shape[0],2*self.ny,2*self.nx))
        for i in range(self.nz):
            expand_v[:,:self.ny,:self.nx] = v_orig[:,i,:,:]
            res += self.xp.fft.fftn(expand_v,axes=[1,2],norm='ortho').astype(np.complex64)*self.eigs[i]
        res = self.xp.real(self.xp.fft.ifft2(res,norm='ortho'))[:,:self.ny,:self.nx]
        res = res.astype(self.xp.float32)
        if v.ndim == 1:
            return res.ravel()
        else:
            return res.reshape(v_orig.shape[0],-1)

    def rmatvec(self,v):
        ''' Apply the transpose of the underlying Toeplitz matrix (i.e. T^T) to
           a vector v.'''
        v_orig = v.reshape(-1,self.ny,self.nx)
        expand_v = np.zeros((v_orig.shape[0],2*self.ny,2*self.nx))
        expand_v[:,:self.ny,:self.nx] = v_orig
        res = self.xp.zeros((v_orig.shape[0],self.nz,self.ny,self.nx))
        v_freq = self.xp.fft.fft2(expand_v,norm='ortho')
        for i in range(self.nz):
            tmp = v_freq * self.xp.conj(self.eigs[i])
            res[:,i,:,:] = self.xp.real(self.xp.fft.ifft2(tmp,norm='ortho'))[:,:self.ny,:self.nx]
        res = res.astype(self.xp.float32)
        if v.ndim == 1:
            return res.ravel()
        else:
            return res.reshape(v_orig.shape[0],-1)

    def cleanup(self):
        self.eigs = None
        self.m_eigs = None
        if self.xp is cupy:
            mempool = cupy.get_default_memory_pool()
            pinned_mempool = cupy.get_default_pinned_memory_pool()
            mempool.free_all_blocks()
            pinned_mempool.free_all_blocks()

class LSQOperator:
    '''This class used to define a Least Square problem by a multilevel
       symmetric Toepliz matrix.'''
    def __init__(self,toep,diag=None,coef=1.,do_precond=False,gtoep=None,*args,**kwargs):
        '''The underlying Toeplitz matrix is determined by it's 1-st row toep.

        Args:
            toep: the 1-st row of the underlying Toeplitz matrix in proper shape.
            diag: the diagonal entries of a diagonal matrix. Used for
                regularization.
            coef: coefficient before the regularization term.
            do_precond: whether or not do preconditioning.
        '''
        if gtoep is None:
            self.gtoep = GToepOperator(toep)
        else:
            self.gtoep = gtoep
        self.xp = self.gtoep.xp
        self.nz,self.ny,self.nx = toep.shape
        self.shape = (np.prod(toep.shape),np.prod(toep.shape))
        self.diag = diag
        self.coef = coef
        self.do_precond = do_precond
        self.dtype = toep.dtype

    def matvec(self,v):
        ''' If do_precond is True, apply T^T*M^{-T}*M^{-1}*T + coef*diag**2 to
            vector v, else apply T^T*T + coef*diag**2 to vector v'''
        tmp = self.gtoep.matvec(v)
        if self.do_precond:
            tmp = self.gtoep.matvec_prec_sym(tmp)
        tmp = self.gtoep.rmatvec(tmp)
        if not self.diag is None:
            tmp += self.coef*self.diag**2*v
        return tmp

    def cleanup(self):
        self.gtoep.cleanup()
        del(self.gtoep)
        self.diag = None
        if self.xp is cupy:
            mempool = cupy.get_default_memory_pool()
            pinned_mempool = cupy.get_default_pinned_memory_pool()
            mempool.free_all_blocks()
            pinned_mempool.free_all_blocks()

def cg(A, b, x=None, tol=1.0e-5, max_iter=None):
    # Note that this function works even tensors 'A' and 'b' are NumPy or CuPy
    # arrays.
    if use_gpu > 0:
        import cupy
        xp = cupy.get_array_module(b)
    else:
        xp = np
    if max_iter is None:
        max_iter = 10*len(b)
    if x is None:
        x = xp.zeros_like(b, dtype=np.float32)
    r0 = b - A.matvec(x)
    p = r0
    now = time.time()
    for i in range(max_iter):
        a = xp.inner(r0, r0) / xp.inner(p, A.matvec(p))
        x += a * p
        r1 = r0 - a * A.matvec(p)
        res = xp.linalg.norm(r1)
        if  res < tol:
            return x
        b = xp.inner(r1, r1) / xp.inner(r0, r0)
        p = r1 + b * p
        r0 = r1
        if time.time() - now > 10:
            now = time.time()
            print('iter: {:8d}  residual: {:16.7e}'.format(i,xp.asnumpy(res)))
            #print('iter: {:8d}  residual: {:16.7e}'.format(i,cupy.asnumpy(res)))
    print('Failed to converge. Increase max-iter or tol.')
    return x
