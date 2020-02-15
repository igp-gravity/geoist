import os
import pathlib
from datetime import datetime
from functools import wraps
from pathos.multiprocessing import Pool
import numpy as np
from scipy import linalg as splin
from scipy import sparse as spsparse
from scipy.optimize import minimize
from scipy.optimize import Bounds
from scipy.optimize import LinearConstraint
import h5py

from geoist import gridder
from geoist.pfm import prism
from geoist.inversion.mesh import PrismMesh
from geoist.inversion import walsh
from geoist.inversion import toeplitz as tptz
from geoist.others import utils

print_level = -1 # control indentation of prints.
last_print_level = -2
use_gpu = 0

# if use_gpu > 0:
#     import cupy as cp
def check_gpu():
    print(use_gpu)    
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
    if use_gpu > 0:
        import cupy as cp
        mempool = cp.get_default_memory_pool()
        pinned_mempool = cp.get_default_pinned_memory_pool()
        mempool.free_all_blocks()
        pinned_mempool.free_all_blocks()
    else:
        print('NO GPU BE USED!!!')


class SmoothOperator:
    def __init__(self):
        self.axis = {'zyx':{'x':-1,'y':-2,'z':-3},
                     'xyz':{'x':-3,'y':-2,'z':-1}}

    def derivation(self,v,component='dx',array_order='zyx'):
        for axis_i in component[1:]:
            slices = [slice(None)]*v.ndim
            slices[self.axis[array_order][axis_i]] = slice(-1,None,-1)
            v = np.diff(v[tuple(slices)],axis=self.axis[array_order][axis_i])[tuple(slices)]
        return v

    def rderivation(self,v,component='dx',array_order='zyx'):
        for axis_i in component[-1:0:-1]:
            slices = [slice(None)]*v.ndim
            slices[self.axis[array_order][axis_i]] = 0
            shape = list(v.shape)
            shape[self.axis[array_order][axis_i]] = 1
            prepend=np.zeros_like(v[tuple(slices)].reshape(tuple(shape)))
            append=np.zeros_like(v[tuple(slices)].reshape(tuple(shape)))
            v = np.diff(v,
                           axis=self.axis[array_order][axis_i],
                           prepend=prepend,
                           append=append)
        return v

class AbicLSQOperator(tptz.LSQOperator):
    '''An operator doing matrix vector multiplication. The matrix is:
        $\alpha_g G^TG + \sum \alpha_i W^TB_i^TB_iW$. Where $\alpha$'s are
        weights, $G$ is kernel matrix, $W$ is depth constraint, $B_i$'s are
        other constrains.
    '''
    def __init__(self,
                 toep,
                 depth_constraint=None,
                 smooth_components=set(),
                 refer_constraint=None,
                 weights=None):
        super().__init__(toep)
        self.weights = weights
        self.depth_constraint = depth_constraint
        self.refer_constraint = refer_constraint
        self.smooth_components = smooth_components
        self.smop = SmoothOperator()
        if self.weights is None:
            self.weights = {'bound':1,'obs':1,'depth':1,'refer':1,'dx':1,'dy':1,'dz':1}

    def matvec(self,v):
        tmp = self.gtoep.matvec(v)
        tmp = self.weights['obs']*self.gtoep.rmatvec(tmp)
        if 'depth' in self.weights.keys():
            v = self.depth_constraint*v
        if 'refer' in self.weights.keys():
            tmp += self.weights['refer']*self.weights['depth']*self.depth_constraint*self.refer_constraint**2*v
        for key in self.smooth_components:
            tmp2 = v.reshape(-1,self.nz,self.ny,self.nx)
            tmp2 = self.smop.derivation(tmp2,component=key)
            tmp2 = self.smop.rderivation(tmp2,component=key)
            if v.ndim == 1:
                tmp += self.weights[key]*self.weights['depth']*self.depth_constraint*tmp2.ravel()
            else:
                tmp += self.weights[key]*self.weights['depth']*self.depth_constraint*tmp2.reshape(v.shape[0],-1)
        return tmp

class AbicLSQOperator2:
    def __init__(self,op_in):
        self.op = op_in

    def matvec(self,v):
        tmp = self.op.matvec(v[:-1]) - self.op.gtoep.rmatvec(np.ones(self.op.nx*self.op.ny))*v[-1]
        tmp2 = -np.sum(self.op.gtoep.matvec(v[:-1]))
        return np.vstack([tmp.reshape(-1,1),np.array(tmp2).reshape(-1,1)])

class GravInvAbicModel:
    def __init__(self,
                 nzyx=[4,4,4],
                 smooth_components=None,
                 depth_constraint=None,
                 model_density=None,
                 refer_density=None,
                 weights=None,
                 source_volume=None,
                 smooth_on='m',
                 subtract_mean=True,
                 data_dir='/data/gravity_inversion'):
        self.gpu_id = 2
        self.subtract_mean = subtract_mean
        self._nz,self._ny,self._nx = nzyx
        self.smooth_on = smooth_on
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
        if smooth_components is None:
            self._smooth_components = (set(weights.keys()) - set(['depth',
                                                                  'obs',
                                                                  'bound',
                                                                  'refer']))
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
        self.smop = SmoothOperator()
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

    def gen_kernel(self, process = 1):
        def calc_kernel(i):
            return prism.gz(self.xp[0:1],self.yp[0:1],self.zp[0:1],[self.mesh[i]])
        if process > 1: #Winodws MP running has possiblely errors.
            print('Number of process:',process)
            with Pool(processes=process) as pool:
                kernel0 = pool.map(calc_kernel,range(len(self.mesh)))
        else:
           kernel0 = [calc_kernel(i) for i in range(len(self.mesh))]

        self.kernel0 = np.array(kernel0).reshape(self.nz,self.ny,self.nx)
        self.kernel_op = AbicLSQOperator(self.kernel0,
                                         depth_constraint=self.constraints['depth'],
                                         smooth_components=self._smooth_components,
                                         refer_constraint=self.constraints['refer'],
                                         weights=self._weights)

    def _diagvec(self,vec=None,diag=None):
        if vec.ndim == 1:
            return vec * diag
        else:
            return  vec * diag.reshape(1,-1)

    @timeit
    def walsh_transform(self,keys=None):
        if keys is None:
            keys = ['kernel'] + list(self.constraints.keys()) + list(self._smooth_components)
        else:
            keys = keys
        
        if use_gpu > 0:
            import cupy as cp
            
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
                        if key == 'kernel':
                            res = f['kernel']['source_volume'][:] - np.array(self.source_volume)
                            res = np.linalg.norm(res)/np.linalg.norm(np.array(self.source_volume))
                            if res > 1.0e-3:
                                is_stored[key] = False
                    except KeyError:
                        continue
        self._gen_walsh_matrix()
        logn = int(np.ceil(np.log2(self._nx*self._ny*self._nz)))
        norm_walsh = 1./(np.sqrt(2)**logn)
        blocks = ['0','1','2','3']
        matvec_op = {'kernel':self.kernel_op.gtoep.matvec,
                     'depth': lambda x: self._diagvec(x,diag=np.sqrt(self.constraints['depth']))
                 }
        for key in self._smooth_components:
            matvec_op[key] = lambda x: self.smop.derivation(x.reshape(-1,self.nz,self.ny,self.nx),component=key).reshape(x.shape[0],-1)
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
                    
                    if use_gpu > 0:
                        with cp.cuda.Device(self.gpu_id):
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
                    else:
                        res = np.zeros((step,step))
                        j = 0
                        while j*step < part_walsh.shape[1]:
                            tmp_block_gpu = np.asarray(part_walsh[:,j*step:(j+1)*step])
                            res += tmp_block_gpu @ tmp_block_gpu.T
                            j += 1
                        if key in self._smooth_components:
                            res[np.abs(res)<1.0e-1*norm_walsh] = 0.     
                            
                    dxyz_group.create_dataset(blocks[i],data=res)
        if ('depth' in keys) and (not is_stored['depth']):
            with h5py.File(self.fname,mode='a') as f:
                try:
                    del f['depth']['constraint']
                except KeyError:
                    pass
                dxyz_group = f['depth']
                dxyz_group.create_dataset('constraint',data=self.constraints['depth'])
        if ('kernel' in keys) and (not is_stored['kernel']):
            with h5py.File(self.fname,mode='a') as f:
                try:
                    del f['kernel']['source_volume']
                except KeyError:
                    pass
                dxyz_group = f['kernel']
                dxyz_group.create_dataset('source_volume',data=np.array(self._source_volume))

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
            for key in self._smooth_components:
                tmp2 = v.reshape(-1,self._nz,self._ny,self._nx)
                tmp2 = self.smop.derivation(tmp2,component=key)
                tmp2 = self.smop.rderivation(tmp2,component=key)
                if v.ndim == 1:
                    self.rhs += self._weights[key]*self._weights['depth']*self.constraints['depth']*tmp2.ravel()
                else:
                    self.rhs += self._weights[key]*self._weights['depth']*self.constraints['depth']*tmp2.reshape(v.shape[0],-1)

    @timeit
    def do_linear_solve(self):
        self.do_linear_solve_quiet()

    def do_linear_solve_quiet(self):
        self._gen_rhs()
        if self.subtract_mean:
            sum_obs = np.sum(self.obs_data)
            tmp_b = np.zeros(len(self.rhs)+1) #np.zeros(len(self.rhs+1)) #chenshi
            tmp_b[:-1] = self.rhs
            tmp_b[-1] = sum_obs
            tmp_op = AbicLSQOperator2(self.kernel_op)
            self.solution = spsparse.linalg.cg(tmp_op,tmp_b,tol=1.0e-5)[0]
        else:
            self.solution = spsparse.linalg.cg(self.kernel_op,self.rhs,tol=1.0e-5)[0]

    @timeit
    def calc_u(self,solved=False,x=None):
        return self.calc_u_quiet(solved,x)

    @timeit
    def calc_min_u(self,solved=False,x=None):
        return self.calc_u_quiet(solved,x)

    def calc_u_quiet(self,solved=False,x=None):
        if x is None:
            if not solved:
                self.do_linear_solve_quiet()
            x = self.solution
        self.min_u_val = self._weights['obs']*np.linalg.norm(self.kernel_op.gtoep.matvec(x) - self.obs_data)**2
        if ('refer' in self._weights.keys()) and (self.smooth_on == 'm-m0'):
            v = x - self.constraints_val['refer']
        else:
            v = x
        if 'depth' in self._weights.keys():
            v = np.sqrt(self._weights['depth'])*self.constraints['depth']*v
        for key in self._smooth_components:
            tmp2 = self.smop.derivation(v.reshape(self._nz,self._ny,self._nx),
                                        component=key)
            self.min_u_val += self._weights[key]*np.linalg.norm(tmp2.ravel())**2
        if 'refer' in self._weights.keys():
            v = x - self.constraints_val['refer']
            if 'depth' in self._weights.keys():
                v = np.sqrt(self._weights['depth'])*self.constraints['depth']*v
            self.min_u_val += self._weights['refer'] *np.linalg.norm(v)**2
        return self.min_u_val

    def jac_u(self,x=None):
        res = self.kernel_op.matvec(x) - self.rhs
        return 2.*res

    def hessp_u(self,x,v):
        res = self.kernel_op.matvec(v)
        return 2.*res

    @timeit
    def bound_optimize(self,x0=None):
        density_bounds = Bounds(self.min_density,self.max_density)
        if x0 is None:
            x0 = np.zeros(self._nx*self._ny*self._nz)+(self.max_density - self.min_density)/2.

        self.bound_solution = minimize(lambda x:self.calc_u_quiet(solved=True,x=x),
                                       x0,
                                       method='trust-constr',
                                       jac=self.jac_u,
                                       hessp=self.hessp_u,
                                       bounds=density_bounds,
                                       )

    def lasso_target(self,x):
#        self.min_u_val = self._weights['obs']*np.linalg.norm(self.kernel_op.gtoep.matvec(x) - self.obs_data)**2
        self.min_u_val = self._weights['obs']*np.linalg.norm(self.kernel_op.gtoep.matvec(x) - self.obs_data)
        if ('refer' in self._weights.keys()) and (self.smooth_on == 'm-m0'):
            v = x - self.constraints_val['refer']
        else:
            v = x
        if 'depth' in self._weights.keys():
            v = np.sqrt(self._weights['depth'])*self.constraints['depth']*v
        for key in self._smooth_components:
            tmp2 = self.smop.derivation(v.reshape(self._nz,self._ny,self._nx),
                                        component=key)
            self.min_u_val += self._weights[key]*np.linalg.norm(tmp2.ravel())
        if 'refer' in self._weights.keys():
            v = x - self.constraints_val['refer']
            if 'depth' in self._weights.keys():
                v = np.sqrt(self._weights['depth'])*self.constraints['depth']*v
            self.min_u_val += self._weights['refer'] *np.linalg.norm(v)
        return self.min_u_val

    def lasso_jac(self,x):
#        jac = self.kernel_op.gtoep.rmatvec(self.kernel_op.gtoep.matvec(x)) - self.kernel_op.gtoep.rmatvec(self.obs_data)
#        jac = 2.0*self._weights['obs']*jac
        jac = self.kernel_op.gtoep.rmatvec(self.kernel_op.gtoep.matvec(x)) - self.kernel_op.gtoep.rmatvec(self.obs_data)
        jac = self._weights['obs']*jac
        norm_res = np.linalg.norm(self.obs_data - self.kernel_op.gtoep.matvec(x))
        jac = jac/norm_res
        if 'refer' in self.weights.keys():
            v = x - self.constraints_val['refer']
            if 'depth' in self._weights.keys():
                v = self.constraints['depth']*v
            norm_refer = np.linalg.norm(v)
            jac += self.weights['refer']*self.weights['depth']*self.constraints['depth']*v/norm_refer

        for key in self.smooth_components:
            v = x
            if 'depth' in self._weights.keys():
                v = self.constraints['depth']*v
            tmp2 = self.smop.derivation(v.reshape(self._nz,self._ny,self._nx),
                                        component=key)
            smooth_norm = np.linalg.norm(tmp2.ravel())
            tmp2 = self.smop.rderivation(tmp2,component=key)
            jac += self.weights[key]*self.weights['depth']*self.constraints['depth']*tmp2.ravel()/smooth_norm
        return jac

    def lasso_hessp(self,x,v):
#        res = self.kernel_op.gtoep.rmatvec(self.kernel_op.gtoep.matvec(v))
#        res = 2.0*self._weights['obs']*res
        norm_res = np.linalg.norm(self.obs_data - self.kernel_op.gtoep.matvec(x))
        res = self.kernel_op.gtoep.rmatvec(self.kernel_op.gtoep.matvec(v))/norm_res
        gradient_res = (self.kernel_op.gtoep.rmatvec(self.kernel_op.gtoep.matvec(x))
                        - self.kernel_op.gtoep.rmatvec(self.obs_data))
        res -= np.dot(gradient_res,v)/norm_res**3 * gradient_res
        res *= self._weights['obs']
        if 'refer' in self.weights.keys():
            v2 = x - self.constraints_val['refer']
            if 'depth' in self._weights.keys():
                v2 = self.constraints['depth']*v2
            norm_refer = np.linalg.norm(v2)
            res += self.weights['refer']*self.weights['depth']/norm_refer*self.constraints['depth']*v
            grad_ref = self.constraints['depth']*v2
            res -= (self.weights['refer']*self.weights['depth']
                    *(np.dot(v,grad_ref)/norm_refer**3*grad_ref))
        for key in self.smooth_components:
            v2 = x
            if 'depth' in self._weights.keys():
                v2 = self.constraints['depth']*v2
            tmp2 = self.smop.derivation(v2.reshape(self._nz,self._ny,self._nx),
                                        component=key)
            smooth_norm = np.linalg.norm(tmp2.ravel())
            tmp2 = self.smop.rderivation(tmp2,component=key)
            grad_sm = self.constraints['depth']*tmp2.ravel()

            tmp2 = v.reshape(self.nz,self.ny,self.nx)
            tmp2 = self.smop.derivation(tmp2,component=key)
            tmp2 = self.smop.rderivation(tmp2,component=key)
            tmp2 = self.constraints['depth']*tmp2.ravel()
            res += (self._weights['depth']*self._weights[key]/smooth_norm
                    *(tmp2 - np.dot(v,grad_sm)*grad_sm/smooth_norm**2))
        return res


    @timeit
    def lasso_optimize(self,x0=None):
        density_bounds = Bounds(self.min_density,self.max_density)
        if x0 is None:
            x0 = (np.random.rand(self._nx*self._ny*self._nz)
                  *(self.max_density - self.min_density)/2.
                  +(self.max_density + self.min_density)/2.)
        self.bound_solution = minimize(self.lasso_target,
                                       x0,
                                       method='trust-constr',
                                       jac=self.lasso_jac,
                                       hessp=self.lasso_hessp,
                                       bounds=density_bounds
                                       )

    def calc_res(self):
        self.residuals = dict()
        self.stds = dict()
        self.residuals['obs'] = np.linalg.norm(self.kernel_op.gtoep.matvec(self.solution)-self.obs_data)**2
        self.stds['obs'] = np.std(self.kernel_op.gtoep.matvec(self.solution)-self.obs_data)
        for key in self._smooth_components:
            tmp2 = self.solution.reshape(self._nz,self._ny,self._nx)
            if ('refer' in self.constraints_val.keys()) and (self.smooth_on == 'm-m0'):
                tmp2 -= self.constraints_val['refer'].reshape(self._nz,self._ny,self._nx)
            tmp2 = self.smop.derivation(tmp2,component=key)
            self.residuals[key] = np.linalg.norm(tmp2.ravel())**2
            self.stds[key] = np.std(tmp2.ravel())
        if 'refer' in self.constraints_val.keys():
            self.residuals['refer'] = np.linalg.norm(self.solution.ravel()-self.constraints_val['refer'].ravel())**2
            self.stds['refer'] = np.std(self.solution.ravel()-self.constraints_val['refer'].ravel())

    def calc_log_prior_total_det_quiet(self):
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
                if use_gpu > 0:
                    import cupy as cp
                    with cp.cuda.Device(self.gpu_id):
                        tmp_block_gpu = cp.asarray(tmp_block,dtype=np.float32)
                        eigs = cp.linalg.eigvalsh(tmp_block_gpu)
                        prior_eigs[i_b*step:(i_b+1)*step] = cp.asnumpy(eigs)
                        self.log_prior_det_val += cp.asnumpy(cp.sum(cp.log(eigs)))
                        tmp_block_gpu = None
                        eigs = None
                        free_gpu()
                    tmp_block += self._weights['obs']*f['kernel'][block][:]
                    with cp.cuda.Device(self.gpu_id):
                        tmp_block_gpu = cp.asarray(tmp_block,dtype=np.float32)
                        eigs = cp.linalg.eigvalsh(tmp_block_gpu)
                        total_eigs[i_b*step:(i_b+1)*step] = cp.asnumpy(eigs)
                        self.log_total_det_val += cp.asnumpy(cp.sum(cp.log(eigs)))
                        tmp_block_gpu = None
                        eigs = None
                        free_gpu()
                else:
                    tmp_block_gpu = np.asarray(tmp_block,dtype=np.float32)
                    eigs = np.linalg.eigvalsh(tmp_block_gpu)
                    prior_eigs[i_b*step:(i_b+1)*step] = eigs
                    self.log_prior_det_val += np.sum(np.log(eigs))
                    tmp_block_gpu = None
                    eigs = None
                    tmp_block += self._weights['obs']*f['kernel'][block][:]
                    tmp_block_gpu = np.asarray(tmp_block,dtype=np.float32)
                    eigs = np.linalg.eigvalsh(tmp_block_gpu)
                    total_eigs[i_b*step:(i_b+1)*step] = eigs
                    self.log_total_det_val += np.sum(np.log(eigs))
                    tmp_block_gpu = None
                    eigs = None                    

                       
        if use_gpu > 0:            
            self.log_prior_det_val = cp.asnumpy(self.log_prior_det_val)
            self.log_total_det_val = cp.asnumpy(self.log_total_det_val)
        else:
            self.log_prior_det_val = self.log_prior_det_val
            self.log_total_det_val = self.log_total_det_val
            
        self.eigs = {'prior':prior_eigs,'total':total_eigs}
        return self.log_prior_det_val,self.log_total_det_val

    @timeit
    def calc_log_prior_total_det(self):
        return self.calc_log_prior_total_det_quiet()

    def calc_log_obs_det_quiet(self):
        self.log_obs_det_val = np.log(self._weights['obs'])*len(self.obs_data)
        return self.log_obs_det_val

    @timeit
    def calc_log_obs_det(self):
        return self.calc_log_obs_det_quiet()

    @timeit
    def calc_abic(self):
        '''-log_prior_det_value+log_total_det-log_obs_det+min_u'''
        self.calc_log_prior_total_det()
        self.calc_u()
        self.calc_log_obs_det()
        self.abic_val = (self.log_total_det_val
                         + self.min_u_val
                         - self.log_prior_det_val
                         - self.log_obs_det_val)
        return self.abic_val

    def calc_abic_quiet(self):
        '''-log_prior_det_value+log_total_det-log_obs_det+min_u'''
        self.calc_log_prior_total_det_quiet()
        self.calc_u_quiet()
        self.calc_log_obs_det_quiet()
        self.abic_val = (self.log_total_det_val
                         + self.min_u_val
                         - self.log_prior_det_val
                         - self.log_obs_det_val)
        return self.abic_val

    def _abic_optimize_exp(self):
        #optimize_keys = list(set(self._weights.keys())-set(['depth']))
        optimize_keys = list(self._weights.keys())
        def abic_target(x):
            for i,key in enumerate(optimize_keys):
                self._weights[key] = np.exp(x[i])
            return self.calc_abic_quiet()
        x0 = np.zeros(len(self._weights))
        for i,key in enumerate(optimize_keys):
            x0[i] = np.log(self._weights[key])
        self.abic_optimize_summary = minimize(abic_target,
                                              x0,
                                              method='Nelder-Mead')

    def _abic_optimize_bound(self):
        optimize_keys = list(self._weights.keys())
        def abic_target(x):
            for i,key in enumerate(optimize_keys):
                self._weights[key] = x[i]
            return self.calc_abic_quiet()
        x0 = np.zeros(len(self._weights))
        for i,key in enumerate(optimize_keys):
            x0[i] = self._weights[key]
        weight_constraint = LinearConstraint(np.eye(len(optimize_keys)),0.,np.inf)
        self.abic_optimize_summary = minimize(abic_target,
                                              x0,
                                              method='COBYLA',
                                              constraints=weight_constraint)

    @timeit
    def abic_optimize(self):
        self._abic_optimize_bound()

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
