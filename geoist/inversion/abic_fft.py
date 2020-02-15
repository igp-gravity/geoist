
import os
import pathlib
from datetime import datetime
from functools import wraps
from pathos.multiprocessing import Pool
import numpy as np
from scipy import linalg as splin
from scipy import sparse as spsparse
from scipy.optimize import minimize
import h5py

import cupy as cp

from geoist import gridder
from geoist.pfm import prism
from geoist.inversion.mesh import PrismMesh
from geoist.others import walsh
from geoist.others import toeplitz as tptz
from geoist.others import utils

print_level = -1 # control indentation of prints.
last_print_level = -2

# A helper decorator print time consumption of f.
def timeit(f):
    @wraps(f)
    def wrap(*args,**kwargs):
        global print_level
        global last_print_level
        print_level += 1
        if print_level == last_print_level:
            print('')
        print(' '*4*print_level+'calling {}'.format(f.__name__))
        st = datetime.now()
        res = f(*args,**kwargs)
        ed = datetime.now()
        print(' '*4*print_level+'{} completed in {}'.format(f.__name__,ed-st))
        last_print_level = print_level
        print_level -= 1
        return res
    return wrap

def free_gpu():
    '''free up gpu memory consumption'''
    mempool = cp.get_default_memory_pool()
    pinned_mempool = cp.get_default_pinned_memory_pool()
    mempool.free_all_blocks()
    pinned_mempool.free_all_blocks()

class SmoothOperator:
    def __init__(self,reverse=False):
        self.axis = {'x':-1,'y':-2,'z':-3}
        if reverse:
            self.axis = {'x':-3,'y':-2,'z':-1}

    def diff(self,v,along='dx'):
        for axis_i in axis_list[1:]:
            slices = [slice(None)]*v.ndim
            slices[self.axis[axis_i]] = slice(-1,None,-1)
            return np.diff(v[tuple(slices)],axis=self.axis[axis_i])

    def rdiff(self,v,along='dx'):
        for axis_i in axis_list[1:]:
            slices = [slice(None)]*v.ndim
            slices[self.axis[axis_i]] = 0
            shape = [-1]*v.ndim
            shape[self.axis[axis_i]] = 1
            prepend=np.zeros_like(v[tuple(slices)].reshape(tuple(shape)))
            append=np.zeros_like(v[tuple(slices)].reshape(tuple(shape)))
            return np.diff(v,
                           axis=self.axis[axis_i],
                           prepend=prepend,
                           append=append)


class AbicLSQOperator(tptz.LSQOperator):
    '''An operator doing matrix vector multiplication. The matrix is:
        $\alpha_g G^TG + \sum \alpha_i W^TB_i^TB_iW$. Where $\alpha$'s are
        weights, $G$ is kernel matrix, $W$ is depth constraint, $B_i$'s are
        other constrains.
    '''
    def __init__(self,
                 toep,
                 depth_constraint=None,
                 dxyz_constraint=None,
                 refer_constraint=None,
                 weights=None):
        super().__init__(toep)
        self.weights = weights
        self.depth_constraint = depth_constraint
        self.refer_constraint = refer_constraint
        self.dxyz_constraint = dxyz_constraint
        if self.weights is None:
            self.weights = {'bound':1,'obs':1,'depth':1,'refer':1,'dx':1,'dy':1,'dz':1}

    def matvec(self,v):
        tmp = self.gtoep.matvec(v)
        tmp = self.weights['obs']*self.gtoep.rmatvec(tmp)
        if 'depth' in self.weights.keys():
            v = self.depth_constraint*v
        if 'refer' in self.weights.keys():
            tmp += self.weights['refer']*self.weights['depth']*self.depth_constraint*self.refer_constraint**2*v
        if not self.dxyz_constraint is None:
            spaces = {'dz':self.nx*self.ny*(self.nz-1),
                      'dy':self.nx*(self.ny-1),
                      'dx':self.nx-1}
            for key,constraint in self.dxyz_constraint.items():
                if not key in self.weights.keys():
                    continue
                tmp2 = v.reshape(-1,*constraint.shape)
                fft_comp = list(range(tmp2.ndim))[1:]
                tmp2 = self.xp.fft.ifftn(self.xp.fft.fftn(tmp2,axes=fft_comp)*constraint,axes=fft_comp).real
                slices = [slice(None)]*tmp2.ndim
                slices[-1] = slice(spaces[key],None)
                tmp2[tuple(slices)] = 0
                tmp2 = self.xp.real(self.xp.fft.ifftn(self.xp.fft.fftn(tmp2,axes=fft_comp)*self.xp.conj(constraint),axes=fft_comp))
                if v.ndim == 1:
                    tmp += self.weights[key]*self.weights['depth']*self.depth_constraint*tmp2.ravel()
                else:
                    tmp += self.weights[key]*self.weights['depth']*self.depth_constraint*tmp2.reshape(v.shape[0],-1)
        return tmp

class GravInvAbicModel:
    def __init__(self,
                 nzyx=[4,4,4],
                 smooth_components=['dx','dy','dz'],
                 depth_constraint=None,
                 model_density=None,
                 refer_density=None,
                 weights=None,
                 source_volume=None,
                 smooth_on='m',
                 data_dir='/data/gravity_inversion'):
        self._nz,self._ny,self._nx = nzyx
        self.smooth_on = smooth_on
        self.dxyz_shapes = {'dx':(self._nz,self._ny,self._nx),
                  'dy':(self._nz,self._nx*self._ny),
                  'dz':(self._nx*self._ny*self._nz,)}
        self.dxyz_spaces = {'dx':self._nx-1,
                  'dy':self._nx*(self._ny-1),
                  'dz':self._nx*self._ny*(self._nz-1)}
        self.data_dir = data_dir
        self.gen_model_name()
        self.nobsx = nzyx[2]
        self.nobsy = nzyx[1]
        self.source_volume = source_volume
        if model_density is None:
            self._model_density = None
        else:
            self._model_density = model_density.ravel()
        self._smooth_components = smooth_components
        self.constraints = dict()
        self.constraints_val = dict()
        if depth_constraint is None:
            self.constraints['depth'] = np.ones(np.prod(nzyx))
            self.constraints_val['depth'] = None
        else:
            self.constraints['depth'] = (depth_constraint.reshape(-1,1)*np.ones((1,self._nx*self._ny))).ravel()
            self.constraints_val['depth'] = 0
        if refer_density is None:
            self.constraints['refer'] = None
            self.constraints_val['refer'] = None
        else:
            self.constraints['refer'] = np.ones(self._nx*self._ny*self._nz)
            self.constraints_val['refer'] = refer_density.ravel()
        self._weights = weights
        if not 'depth' in self._weights.keys():
            self._weights['depth'] = 1.0
        self._gen_dxyz_constraint()
        self.kernel_op = None
        self.abic_val = 0
        self.log_total_det_val = 0
        self.log_prior_det_val = 0
        self.log_obs_det_val = 0
        self.min_u_val = 0
        self.min_density = -1.0e4
        self.max_density = 1.0e4

    @property
    def source_volume(self):
        return self._source_volume
    @source_volume.setter
    def source_volume(self,value):
        self._source_volume = value
        self.gen_mesh()

    def gen_model_name(self):
        self.model_name = '{}x{}x{}'.format(self._nx,self._ny,self._nz)
        self.fname = pathlib.Path(self.data_dir)/pathlib.Path(self.model_name+'.h5')

    @property
    def weights(self):
        return self._weights
    @weights.setter
    def weights(self,values):
        self._weights = values
        if not self.kernel_op is None:
            self.kernel_op.weights = self._weights

    @property
    def smooth_components(self):
        return self._smooth_components
    @smooth_components.setter
    def smooth_components(self,values):
        self._smooth_components = values
        self._gen_dxyz_constraint()
        if not self.kernel_op is None:
            self.kernel_op.dxyz_constraint = self.dxyz_constraint

    @timeit
    def _gen_dxyz_constraint(self):
        '''first generate multi-level circulant matrix, constraint of dx is a part of it. then calculate it's eigenvalues.
        self._dx stores the eigenvalues finally. When multiply it with a vector, specific element should be discarded'''
        self.dxyz_constraint = dict()
        for component in self._smooth_components:
            tmp = np.zeros(self.nx*self.ny*self.nz)
            tmp[0] = 1
            tmp[self.dxyz_spaces[component]] = -1
            tmp = tmp.reshape(self.dxyz_shapes[component])
            self.dxyz_constraint[component] = np.fft.fftn(tmp)
            self.constraints[component] = self.dxyz_constraint[component]

    @property
    def refer_density(self):
        return self.constraints_val['refer'].reshape(self._nz,self._ny,self._nx)
    @refer_density.setter
    def refer_density(self,value):
        self.constraints_val['refer'] = value.ravel()

    @property
    def nx(self):
        return self._nx
    @nx.setter
    def nx(self,value):
        self._nx = value
        self.nobsx = self._nx
        self.gen_model_name()
        if not self.constraints['depth'] is None:
            self.constraints['depth'] = self.constraints['depth'].reshape(self._nz,-1)[:,0]*np.ones((1,self._nx*self._ny))
            self.constraints['depth'] = self.constraints['depth'].ravel()
        self.constraints['refer'] = np.ones(self._nx*self._ny*self._nz)

    @property
    def ny(self):
        return self._ny
    @ny.setter
    def ny(self,value):
        self._ny = value
        self.nobsy = self._ny
        self.gen_model_name()
        if not self.constraints['depth'] is None:
            self.constraints['depth'] = self.constraints['depth'].reshape(self._nz,-1)[:,0]*np.ones((1,self._nx*self._ny))
            self.constraints['depth'] = self.constraints['depth'].ravel()
        self.constraints['refer'] = np.ones(self._nx*self._ny*self._nz)

    @property
    def nz(self):
        return self._nz
    @nz.setter
    def nz(self,value):
        self._nz = value
        self.gen_model_name()
        self.constraints['refer'] = np.ones(self._nx*self._ny*self._nz)
        print("Warning: nz changed. \nDon't forget setting depth constraints.")

    @property
    def model_density(self):
        return(self._model_density.reshape(self.nz,self.ny,self.nx))
    @model_density.setter
    def model_density(self,value):
        self._model_density = value.ravel()

    def gen_mesh(self,height = -1):
        shape = (self._nz, self._ny, self._nx)
        self.mesh = PrismMesh(self._source_volume, shape)
        density = np.ones(shape)*1.0e3
        self.mesh.addprop('density', density.ravel())
        # generate obs grid
        # coordinate: x North-South,y East-West
        # gridder is in the order: (nx,ny)
        self.obs_area = (self._source_volume[0]+0.5*self.mesh.dims[0],
                         self._source_volume[1]-0.5*self.mesh.dims[0],
                         self._source_volume[2]+0.5*self.mesh.dims[1],
                         self._source_volume[3]-0.5*self.mesh.dims[1])
        obs_shape = (self.nobsx, self.nobsy)
        self.xp, self.yp, self.zp = gridder.regular(self.obs_area, obs_shape, z=height)

    def _gen_walsh_matrix(self):
        print('generating walsh_matrix')
        if os.path.exists(self.fname):
            with h5py.File(self.fname,mode='r') as f:
                if not 'walsh_matrix' in f.keys():
                    have_walsh_matrix = False
                else:
                    have_walsh_matrix = True
        else:
            have_walsh_matrix = False
        if have_walsh_matrix:
            return
        walsh_matrix = walsh.walsh_matrix(self._nx*self._ny*self._nz,
                                          normalized=True,
                                          ordering='sequence2',
                                          nxyz=(self._nx,self._ny,self._nz))
        walsh_matrix = walsh_matrix.astype(np.float32)
        step = self._nx*self._ny*self._nz//4
        components = ['0','1','2','3']
        with h5py.File(self.fname,mode='a') as f:
            fgroup = f.create_group('walsh_matrix')
            for i in range(4):
                fgroup.create_dataset(components[i],data=walsh_matrix[i*step:(i+1)*step,:])

    def gen_kernel(self):
        def calc_kernel(i):
            return prism.gz(self.xp[0:1],self.yp[0:1],self.zp[0:1],[self.mesh[i]])
        with Pool(processes=16) as pool:
            kernel0 = pool.map(calc_kernel,range(len(self.mesh)))
        self.kernel0 = np.array(kernel0).reshape(self.nz,self.ny,self.nx)
        self.kernel_op = AbicLSQOperator(self.kernel0,
                                         depth_constraint=self.constraints['depth'],
                                         dxyz_constraint=self.dxyz_constraint,
                                         refer_constraint=self.constraints['refer'],
                                         weights=self._weights)

    def _dxyzvec(self,vec=None,key=None):
        res = vec.reshape(-1,*self.dxyz_shapes[key])
        axes = np.arange(1,res.ndim)
        res = np.fft.ifftn(np.fft.fftn(res,axes=axes)*self.dxyz_constraint[key],axes=axes).real
        slices = [slice(None)]*res.ndim
        slices[-1] = slice(0,self.dxyz_spaces[key])
        if vec.ndim == 1:
            return res[tuple(slices)].ravel()
        else:
            return res[tuple(slices)].reshape(vec.shape[0],-1)

    def _diagvec(self,vec=None,diag=None):
        if vec.ndim == 1:
            return vec * diag
        else:
            return  vec * diag.reshape(1,-1)

    @timeit
    def walsh_transform(self,keys=None):
        if keys is None:
            keys = ['kernel'] + list(self.constraints.keys())
        else:
            keys = keys
        is_stored = dict()
        for key in keys:
            is_stored[key] = False
        if os.path.exists(self.fname):
            with h5py.File(self.fname,mode='r') as f:
                for key in keys:
                    try:
                        if '3' in f[key].keys():
                            is_stored[key] = True
                        if key == 'depth':
                            res = f['depth']['constraint'][:] - self.constraints['depth']
                            res = np.linalg.norm(res)/np.linalg.norm(self.constraints['depth'])
                            if res > 1.0e-3:
                                is_stored[key] = False
                    except KeyError:
                        continue
        self._gen_walsh_matrix()
        logn = int(np.ceil(np.log2(self._nx*self._ny*self._nz)))
        norm_walsh = 1./(np.sqrt(2)**logn)
        blocks = ['0','1','2','3']
        matvec_op = {'kernel':self.kernel_op.gtoep.matvec,
                  'dx': lambda x: self._dxyzvec(x,key='dx'),
                  'dy': lambda x: self._dxyzvec(x,key='dy'),
                  'dz': lambda x: self._dxyzvec(x,key='dz'),
                  'refer': lambda x: self._diagvec(x,diag=self.constraints['refer']),
                  'depth': lambda x: self._diagvec(x,diag=np.sqrt(self.constraints['depth']))
                 }
        is_stored['refer'] = True
        for key in keys:
            if is_stored[key]:
                print('walsh transformation of {} already exists.'.format(key))
                continue
            print('performing walsh transformation on {}.'.format(key))
            step = self.nx*self.ny*self.nz // 4
            if key == 'depth':
                step = self._nz
            with h5py.File(self.fname,mode='a') as f:
                try:
                    del f[key]
                except KeyError:
                    pass
                dxyz_group = f.create_group(key)
                walsh_group = f['walsh_matrix']
                for i in range(4):
                    print("\t progress {}/4".format(i))
                    part_walsh = walsh_group[blocks[i]][:]
                    if key == 'depth':
                        part_walsh = walsh_group[blocks[i]][:self._nz]
                    part_walsh = matvec_op[key](part_walsh)
                    with cp.cuda.Device(2):
                        res = cp.zeros((step,step))
                        j = 0
                        while j*step < part_walsh.shape[1]:
                            tmp_block_gpu = cp.asarray(part_walsh[:,j*step:(j+1)*step])
                            res += tmp_block_gpu @ tmp_block_gpu.T
                            j += 1
                        res = cp.asnumpy(res)
                        if key in self._smooth_components:
                            res[np.abs(res)<1.0e-1*norm_walsh] = 0.
                        tmp_block_gpu = None
                        mempool = cp.get_default_memory_pool()
                        pinned_mempool = cp.get_default_pinned_memory_pool()
                        mempool.free_all_blocks()
                        pinned_mempool.free_all_blocks()
                    dxyz_group.create_dataset(blocks[i],data=res)
        if ('depth' in keys) and (not is_stored['depth']):
            with h5py.File(self.fname,mode='a') as f:
                try:
                    del f['depth_constraint']
                except KeyError:
                    pass
                dxyz_group = f['depth']
                dxyz_group.create_dataset('constraint',data=self.constraints['depth'])

    @property
    def depth_constraint(self):
        return(self.constraints['depth'].reshape(self._nz,-1)[:,0])
    @depth_constraint.setter
    def depth_constraint(self,value):
        self.constraints['depth'] = (value.reshape(-1,1)*np.ones((1,self._nx*self._ny))).ravel()

    @timeit
    def forward(self,model_density=None):
        if model_density is None:
            model_density = self._model_density
        else:
            model_density = model_density.ravel()
        self.obs_data = self.kernel_op.gtoep.matvec(model_density)

    def _gen_rhs(self):
        self.rhs = self._weights['obs']*self.kernel_op.gtoep.rmatvec(self.obs_data)
        if 'depth' in self._weights.keys():
            v = self.constraints['depth']*self.constraints_val['refer']
        if 'refer' in self._weights.keys():
            self.rhs += (self._weights['refer']
                         *self._weights['depth']
                         *self.constraints['depth']
                         *v)
        if self.smooth_on == 'm-m0':
            if not self.dxyz_constraint is None:
                for key,constraint in self.dxyz_constraint.items():
                    if not key in self._weights.keys():
                        continue
                    tmp2 = v.reshape(-1,*constraint.shape)
                    fft_comp = list(range(tmp2.ndim))[1:]
                    tmp2 = np.fft.ifftn(np.fft.fftn(tmp2,axes=fft_comp)*constraint,axes=fft_comp).real
                    slices = [slice(None)]*tmp2.ndim
                    slices[-1] = slice(self.dxyz_spaces[key],None)
                    tmp2[tuple(slices)] = 0
                    tmp2 = np.real(np.fft.ifftn(np.fft.fftn(tmp2,axes=fft_comp)*np.conj(constraint),axes=fft_comp))
                    if v.ndim == 1:
                        self.rhs += self._weights[key]*self._weights['depth']*self.constraints['depth']*tmp2.ravel()
                    else:
                        self.rhs += self._weights[key]*self._weights['depth']*self.constraints['depth']*tmp2.reshape(v.shape[0],-1)

    @timeit
    def do_linear_solve(self):
        self._gen_rhs()
        self.solution = spsparse.linalg.cg(self.kernel_op,self.rhs,tol=1.0e-5)[0]

    @timeit
    def calc_min_u(self,solved=False,x=None):
        if x is None:
            if not solved:
                self.do_linear_solve()
            x = self.solution
        self.min_u_val = self._weights['obs']*np.linalg.norm(self.kernel_op.gtoep.matvec(x) - self.obs_data)**2
        if ('refer' in self._weights.keys()) and (self.smooth_on == 'm-m0'):
            v = x - self.constraints_val['refer']
        else:
            v = x
        if 'depth' in self._weights.keys():
            v = np.sqrt(self._weights['depth'])*self.constraints['depth']*v
        if not self.dxyz_constraint is None:
            for key,constraint in self.dxyz_constraint.items():
                if not key in self._weights.keys():
                    continue
                tmp2 = np.fft.ifftn(
                         np.fft.fftn(v.reshape(constraint.shape))*constraint
                       ).real
                slices = [slice(None)]*constraint.ndim
                slices[-1] = slice(0,self.dxyz_spaces[key])
                self.min_u_val += self._weights[key]*np.linalg.norm(tmp2[tuple(slices)].ravel())**2
        if 'refer' in self._weights.keys():
            v = x - self.constraints_val['refer']
            if 'depth' in self._weights.keys():
                v = np.sqrt(self._weights['depth'])*self.constraints['depth']*v
            self.min_u_val += self._weights['refer'] *np.linalg.norm(v)**2
        return self.min_u_val

    def bound_constraint_u(self,x=None):
        self.calc_min_u(x=x,solved=True)
        log_barrier = np.sum(np.log(x-self.min_density) + np.log(self.max_density-x))
        return self.min_u_val - 2.*self._weights['bound']*log_barrier

    def bound_jac_u(self,x=None):
        res = 0.
        res += self._weights['obs']*(self.kernel_op.gtoep.matvec(x) - self.obs_data)
        if ('refer' in self._weights.keys()) and (self.smooth_on == 'm-m0'):
            v = x - self.constraints_val['refer']
        else:
            v = x
        if 'depth' in self._weights.keys():
            v = self._weights['depth']*self.constraints['depth']*v
        if not self.dxyz_constraint is None:
            for key,constraint in self.dxyz_constraint.items():
                if not key in self._weights.keys():
                    continue
                tmp2 = np.fft.ifftn(
                         np.fft.fftn(v.reshape(constraint.shape))*constraint
                       ).real
                slices = [slice(None)]*constraint.ndim
                slices[-1] = slice(0,self.dxyz_spaces[key])
                res += self._weights[key]*tmp2[tuple(slices)].ravel()
        if 'refer' in self._weights.keys():
            v = x - self.constraints_val['refer']
            if 'depth' in self._weights.keys():
                v = self._weights['depth']*self.constraints['depth']*v
            res += self._weights['refer'] *v
        res += self._weights['bound']*(1./(self.max_density-x) - 1./(x-self.min_density))
        return 2.*res

    def bound_hessp_u(self,x,v):
        res = self.kernel_op.matvec(v)
        hess_diag = 1./(self.max_density-x)**2 + 1./(x-self.min_density)**2
        res += self._weights['bound']*hess_diag*v
        return 2.*res

    def bound_optimize(self,x0=None):
        if x0 is None:
            if 'refer' in self._weights.keys():
                x0 = self.constraints_val['refer']
            else:
                x0 = np.zeros(self._nx*self._ny*self._nz)
        self.solution = minimize(self.bound_constraint_u,
                       x0,
                       method='Newton-CG',
                       jac=self.bound_jac_u,
                       hessp=self.bound_hessp_u)

    def calc_res(self):
        self.residuals = dict()
        self.stds = dict()
        self.residuals['obs'] = np.linalg.norm(self.kernel_op.gtoep.matvec(self.solution)-self.obs_data)**2
        self.stds['obs'] = np.std(self.kernel_op.gtoep.matvec(self.solution)-self.obs_data)
        for key in self.dxyz_constraint.keys():
            try:
                tmp2 = self.solution.reshape(self.dxyz_constraint[key].shape)
                if ('refer' in self.constraints_val.keys()) and (self.smooth_on == 'm-m0'):
                    tmp2 -= self.constraints_val['refer'].reshape(self.dxyz_constraint[key].shape)
                tmp2 = np.fft.ifftn(np.fft.fftn(tmp2)*self.dxyz_constraint[key]).real
                slices = [slice(None)]*tmp2.ndim
                slices[-1] = slice(0,self.dxyz_spaces[key])
                self.residuals[key] = np.linalg.norm(tmp2[tuple(slices)].ravel())**2
                self.stds[key] = np.std(tmp2[tuple(slices)].ravel())
            except KeyError:
                pass
        if 'refer' in self.constraints_val.keys():
            self.residuals['refer'] = np.linalg.norm(self.solution.ravel()-self.constraints_val['refer'].ravel())**2
            self.stds['refer'] = np.std(self.solution.ravel()-self.constraints_val['refer'].ravel())

    @timeit
    def calc_log_prior_total_det(self):
        self.log_prior_det_val = 0
        self.log_total_det_val = 0
        blocks = ['0','1','2','3']
        prior_eigs = np.zeros(self._nx*self._ny*self._nz)
        total_eigs = np.zeros(self._nx*self._ny*self._nz)
        step = self._nx*self._ny*self._nz//4
        try:
            depth_weight = self._weights['depth']
        except KeyError:
            depth_weight = 1.
        with h5py.File(self.fname,mode='r') as f:
            if 'depth' in self._weights.keys():
                depth_walsh = f['depth']['0'][:]
            for i_b,block in enumerate(blocks):
                tmp_block = np.zeros((step,step))
                for dxyz_name in self._smooth_components:
                    try:
                        dxyz_walsh = f[dxyz_name][block][:].reshape(step//self._nz,
                                                                    self._nz,
                                                                    step//self._nz,
                                                                    self._nz)
                        ein_path = np.einsum_path('mi,xiyj,jn->xmyn',
                                                  depth_walsh.T,
                                                  dxyz_walsh,
                                                  depth_walsh,
                                                  optimize='optimal')[0]
                        tmp_multi = np.einsum('mi,xiyj,jn->xmyn',
                                              depth_walsh.T,
                                              dxyz_walsh,
                                              depth_walsh,
                                              optimize=ein_path)
                        tmp_block += depth_weight*self._weights[dxyz_name]*tmp_multi.reshape(step,step)
                    except KeyError:
                        pass
                if 'refer' in self._weights.keys():
                    tmp_multi_small = depth_walsh.T@depth_walsh
                    for i in range(step//self._nz):
                        tmp_block[i*self._nz:(i+1)*self._nz,
                                  i*self._nz:(i+1)*self._nz] += depth_weight*self._weights['refer']*tmp_multi_small
                with cp.cuda.Device(2):
                    tmp_block_gpu = cp.asarray(tmp_block,dtype=np.float32)
                    eigs = cp.linalg.eigvalsh(tmp_block_gpu)
                    prior_eigs[i_b*step:(i_b+1)*step] = cp.asnumpy(eigs)
                    self.log_prior_det_val += cp.asnumpy(cp.sum(cp.log(eigs)))
                    tmp_block_gpu = None
                    eigs = None
                    free_gpu()
                tmp_block += self._weights['obs']*f['kernel'][block][:]
                with cp.cuda.Device(2):
                    tmp_block_gpu = cp.asarray(tmp_block,dtype=np.float32)
                    eigs = cp.linalg.eigvalsh(tmp_block_gpu)
                    total_eigs[i_b*step:(i_b+1)*step] = cp.asnumpy(eigs)
                    self.log_total_det_val += cp.asnumpy(cp.sum(cp.log(eigs)))
                    tmp_block_gpu = None
                    eigs = None
                    free_gpu()
        self.log_prior_det_val = cp.asnumpy(self.log_prior_det_val)
        self.log_total_det_val = cp.asnumpy(self.log_total_det_val)
        self.eigs = {'prior':prior_eigs,'total':total_eigs}
        return self.log_prior_det_val,self.log_total_det_val

    @timeit
    def calc_log_obs_det(self):
        self.log_obs_det_val = np.log(self._weights['obs'])*len(self.obs_data)
        return self.log_obs_det_val

    @timeit
    def calc_abic(self):
        '''-log_prior_det_value+log_total_det-log_obs_det+min_u'''
        self.calc_log_prior_total_det()
        self.calc_min_u()
        self.calc_log_obs_det()
        self.abic_val = (self.log_total_det_val
                         + self.min_u_val
                         - self.log_prior_det_val
                         - self.log_obs_det_val)
        return self.abic_val

    @timeit
    def para_grad(self,x):
        pass

    def u_bound(self):
        pass

    def print_summary(self):
        print('abic values:{}'.format(self.abic_val))
        print('log total det:{}'.format(self.log_total_det_val))
        print('log prior det:{}'.format(self.log_prior_det_val))
        print('log obs det:{}'.format(self.log_obs_det_val))
        print('min u:{}'.format(self.min_u_val))
        print('std:',end=' ')
        print(self.stds)
        print('1/var:',end=' ')
        print({k:1./v**2 for k,v in self.stds.items()})
        print('norms:',end=' ')
        print(self.residuals)
