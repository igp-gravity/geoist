import numpy as np
from geoist.inversion import abic
from geoist import gridder
from geoist.inversion.mesh import PrismMesh
from geoist.pfm import pftrans
from geoist.pfm.giconstants import G, SI2MGAL
import json

def kron_matvec(t, V):
    '''matrix M multiply vectors V where M is kronecker product of A and B,
    i.e. M = A\otimes B.
    Args:
        t (list of ndarray): M = tds[0] \otimes t[1] \otimes t[2] ...
        V (ndarray): vectors to be multiplied.
    Returns:
        res (ndarray): results.
    '''
    shapes = [m.shape[1] for m in t]
    if V.ndim == 1:
        tmp = V.reshape(1,*shapes)
    else:
        tmp = V.reshape(V.shape[0],*shapes)
    n = len(t)
    params = []
    for i,m in enumerate(t):
        params.append(m)
        params.append([i,n+i])
    params.append(tmp)
    params.append([2*n]+list(range(n,2*n)))
    params.append([2*n]+list(range(n)))
    path = np.einsum_path(*params,optimize='optimal')
    res = np.einsum(*params,optimize=path[0])
    res = res.reshape(res.shape[0],-1)
    return res.squeeze()

class FreqInvModel:
    def __init__(self,conf_file=None,**kwargs):
        self.confs = {'nzyx':[4,4,4],
                      'nobsyx':None,
                      'smooth_components':None,
                      'depth_scaling':None,
                      'model_density':None,
                      'refer_densities':None,
                      'weights':None,
                      'source_volume':None,
                      'obs_area':None,
                      'data_dir':'./'}
        confs = dict()
        if not conf_file is None:
            with open(conf_file) as f:
                confs = json.load(f)
        self.confs = {**self.confs,**confs,**kwargs}
        self.nz,self.ny,self.nx = self.confs['nzyx']
        if self.confs['nobsyx'] is None:
            self.nobsx = self.nx
            self.nobsy = self.ny
        else:
            self.nobsy,self.nobsx = self.confs['nobsyx']
        self.source_volume = self.confs['source_volume']
        self.obs_area= self.confs['obs_area']
        if self.confs['model_density'] is None:
            self._model_density = None
        else:
            self._model_density = self.confs['model_density'].ravel()
        self._smooth_components = self.confs['smooth_components']
        if self.confs['depth_scaling'] is None:
            self.depth_scaling = np.ones(self.nz)
        else:
            self.depth_scaling = self.confs['depth_scaling']
        self.refer_densities = self.confs['refer_densities']
        self._weights = self.confs['weights']
        self.smop = abic.SmoothOperator()
        self.kernel_op = None

    def gen_mesh(self,height = 0):
        ''' Generate mesh of the model.
        Args:
            height (float): height of the observations.
        '''
        shape = (self.nz, self.ny, self.nx)
        self.mesh = PrismMesh(self.source_volume, shape)
        self.mesh.addprop('density', self._model_density)
        # generate obs grid
        # coordinate: x North-South,y East-West
        # gridder is in the order: (nx,ny)
        self.gen_obs_grid(height=height)

    def gen_obs_grid(self,height=0):
        if self.obs_area is None:
            self.obs_area = (self.source_volume[0]+0.5*self.mesh.dims[0],
                              self.source_volume[1]-0.5*self.mesh.dims[0],
                              self.source_volume[2]+0.5*self.mesh.dims[1],
                              self.source_volume[3]-0.5*self.mesh.dims[1])
        obs_shape = (self.nobsx, self.nobsy)
        self.xp, self.yp, self.zp = gridder.regular(self.obs_area, obs_shape, z=height)

    def _pad_data(self, data, shape):
        n0 = pftrans._nextpow2(2*shape[0])
        n1 = pftrans._nextpow2(2*shape[1])
        nx, ny = shape
        padx = (n0 - nx)//2
        pady = (n1 - ny)//2
        padded = np.pad(data.reshape(shape), ((padx, padx), (pady, pady)),
                       mode='edge')
        return padded, padx, pady

    def gen_kernel(self,gtype='z'):
        xs = np.array(self.mesh.get_xs())
        ys = np.array(self.mesh.get_ys())
        zs = np.array(self.mesh.get_zs())
        x0 = (xs[:-1] + xs[1:])/2.0
        y0 = (ys[:-1] + ys[1:])/2.0
        a = np.abs(xs[:-1]-xs[1:])/2.0
        b = np.abs(ys[:-1]-ys[1:])/2.0
        nx, ny = self.nobsx, self.nobsy
        xmin, xmax, ymin, ymax = self.obs_area
        dx = (xmax - xmin)/(nx - 1)
        dy = (ymax - ymin)/(ny - 1)

        # Pad the array with the edge values to avoid instability
        shape = (nx,ny)
        self.padded, self.padx, self.pady = self._pad_data(self.zp, shape)
        nxe, nye = self.padded.shape
        M_left=(nxe-nx)/2+1
        M_right=M_left+nx-1
        N_down=(nye-ny)/2+1
        N_up=N_down+ny-1

        XXmin=xmin-dx*(M_left-1)
        XXmax=xmax+dx*(nxe-M_right)
        YYmin=ymin-dy*(N_down-1)
        YYmax=ymax+dy*(nye-N_up)
        # we store kx and ky as 1d array
        self.kx = 2*np.pi*np.array(np.fft.fftfreq(self.padded.shape[0], dx))
        self.ky = 2*np.pi*np.array(np.fft.fftfreq(self.padded.shape[1], dy))
        # kz is 2d array
        self.kz = np.sqrt(np.add.outer(self.ky**2,self.kx**2)).ravel()
        self.kxm = np.ma.array(self.kx, mask= self.kx==0)
        self.kym = np.ma.array(self.ky, mask= self.ky==0)
        self.kzm = np.ma.array(self.kz, mask= self.kz==0)

        complex1 = 0+1j
        ## zs should be depth or coordinate?
        self.C = -8*np.pi*G*SI2MGAL
        self.W = np.exp(self.kz*self.zp[0])
        if gtype == 'zz':
            self.W = self.kz * self.W
        self.dW =((np.exp(-np.outer(self.kz,zs[1:]))
                  -np.exp(-np.outer(self.kz,zs[:-1])))/self.kzm.reshape(-1,1))
        self.dW[self.kzm.mask,:] = 0.
        self.dW[self.kzm.mask,:] += zs[:-1] - zs[1:]
        self.dW = self.dW.data
        self.WX = np.exp(complex1*(XXmin-XXmax+xmax+xmin)*self.kx/2.)/dx
        #self.WX = np.ones_like(self.kx)/dx
        if gtype == 'zx':
            self.WX = complex1 * self.kx
        self.FX = np.sin(np.outer(self.kx,a))/self.kxm.reshape(-1,1)
        self.FX[self.kxm.mask,:] = 0.
        self.FX[self.kxm.mask,:] += a
        self.FX = self.FX * np.exp(-complex1*np.outer(self.kx,x0)) * self.WX.reshape(-1,1)
        self.FX = self.FX.data
        self.WY = np.exp(complex1*(YYmin-YYmax+ymax+ymin)*self.ky/2.)/dy
#        self.WY = np.ones_like(self.ky)/dy
        if gtype == 'zy':
            self.WY = complex1 * self.ky
        self.FY = np.sin(np.outer(self.ky,b))/self.kym.reshape(-1,1)
        self.FY[self.kym.mask,:] = 0.
        self.FY[self.kym.mask,:] += b
        self.FY = self.FY * np.exp(-complex1*np.outer(self.ky,y0)) * self.WY.reshape(-1,1)
        self.FY = self.FY.data
        self.solve_prepared = False

    def _gen_rhs(self):
        # default pad 0s
        self.rhs = self.C*self.W*self.obs_freq.T.ravel()
        self.rhs = self.dW.T*self.rhs
#        self.rhs = kron_matvec([np.eye(self.nz),
#                                np.conj(self.FY.T),
#                                np.conj(self.FX.T)],
#                               self.rhs.ravel())
        self.rhs = kron_matvec([np.conj(self.FY.T),
                                np.conj(self.FX.T)],
                               self.rhs.reshape(self.nz,-1)).ravel()

    def forward(self,v=None,update=False):
        ''' Forward modelling gravity and its frequency representation
        Args:
            v (ndarray): model density, if not given, use self.model_density
            update (bool): If False, self.freq and self.obs_field won't be touched.
                           If True, self.freq and self.obs_field will be updated.
        Returns:
            obs_field (ndarray): gravity.
            freq (ndarray): gravity in frequency domain.
        '''
        if v is None:
            v = self._model_density
        freq = kron_matvec([self.FY,self.FX],v.reshape(self.nz,-1))
        freq = freq * self.dW.T
        freq = self.W * np.sum(freq,axis=0)
        freq = self.C * freq.reshape(self.padded.shape[1],self.padded.shape[0]).T
        obs_field = np.real(np.fft.ifft2(freq))
        obs_field = obs_field[self.padx: self.padx + self.nobsx, self.pady: self.pady + self.nobsy].ravel()
        if update:
            self.freq = freq
            self.obs_field = obs_field
        return freq,obs_field

    def _prepare_solve(self):
        if self.solve_prepared:
            return
        tmpX = self.FX @ np.conj(self.FX.T)
        tmpY = self.FY@np.conj(self.FY.T)
        tmp = np.kron(tmpY,tmpX)
        self.woodSUB = np.zeros_like(tmp)
        for iz in range(self.nz):
            self.woodSUB += self.dW[:,iz].reshape(-1,1)*tmp*np.conj(self.dW[:,iz].reshape(1,-1))
        self.solve_prepared = True

    def do_linear_solve_quiet(self):
        self._gen_rhs()
        if not self.solve_prepared:
            self._prepare_solve()
        weight = np.sum(self._weights['refers'])
        invSUB = self.C**2*weight*self.woodSUB
        invSUB_diag = np.einsum('ii->i',invSUB)
        invSUB_diag += weight**2 * 1./self.W**2
        invSUB = np.linalg.inv(invSUB)
        self.solution = kron_matvec([self.FY,self.FX],self.rhs.reshape(self.nz,-1))
        self.solution = np.sum(self.solution * self.dW.T, axis=0)
        self.solution = invSUB @ self.solution
        self.solution = self.dW.T*self.solution.reshape(1,-1)
        self.solution = kron_matvec([np.conj(self.FY.T),np.conj(self.FX.T)],
                                    self.solution.reshape(self.nz,-1)).ravel()
        self.solution = self.rhs/weight - self.C**2*self.solution

    def calc_res(self):
        pass
