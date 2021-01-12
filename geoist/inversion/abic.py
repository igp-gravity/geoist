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

    def shapes(self,component='dx',model_size=None):
        '''shape information
        Returns:
            shapes(list): shapes[0] is the correct shape of a vector for rderivation operating on.
                          shapes[1] is the shape of derivation matrix.
        '''
        testv = np.zeros(model_size)
        resv = self.derivation(testv,component=component)
        return [resv.shape,(len(resv.ravel()),len(testv.ravel()))]


class AbicLSQOperator(tptz.LSQOperator):
    '''An operator doing matrix vector multiplication. The matrix is:
        :math:`\alpha_g G^TG + \sum \alpha_i W^TB_i^TB_iW`. Where :math:`\alpha`'s are
        weights, :math:`G` is kernel matrix, :math:`W` is depth constraint, :math:`B_i`'s are
        other constrains.
    '''
    def __init__(self,
                 toep,
                 depth_constraint=None,
                 smooth_components=set(),
                 refer_constraints=None,
                 weights=None):
        super().__init__(toep)
        self.weights = weights
        self.depth_constraint = depth_constraint
        self.refer_constraints = refer_constraints
        self.smooth_components = smooth_components
        self.smop = SmoothOperator()
        if self.weights is None:
            self.weights = {'bound':1,'obs':1,'depth':1,'refers':[1],'dx':1,'dy':1,'dz':1}

    def matvec(self,v):
        tmp = self.gtoep.matvec(v)
        tmp = self.weights['obs']*self.gtoep.rmatvec(tmp)
        if 'depth' in self.weights.keys():
            v = self.depth_constraint*v
        if 'refers' in self.weights.keys():
            for refer_weight in self.weights['refers']:
                tmp += refer_weight*self.depth_constraint*v
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
    def __init__(self,conf_file=None,**kwargs):
        self.confs = {'nzyx':[4,4,4],
                      'gpu_id':2,
                      'smooth_components':None,
                      'depth_constraint':None,
                      'model_density':None,
                      'refer_densities':None,
                      'weights':None,
                      'source_volume':None,
                      'smooth_on':'m',
                      'subtract_mean':False,
                      'optimize_keys':None,
                      'mode':'walsh',
                      'data_dir':'/data/gravity_inversion'}
        confs = dict()
        if not conf_file is None:
            with open(conf_file) as f:
                confs = json.load(f)
        self.confs = {**self.confs,**confs,**kwargs}
        self._nz,self._ny,self._nx = self.confs['nzyx']
        self.nobsx = self._nx
        self.nobsy = self._ny
        self.gen_model_name()
        self.source_volume = self.confs['source_volume']
        if self.confs['model_density'] is None:
            self._model_density = None
        else:
            self._model_density = self.confs['model_density'].ravel()
        self._smooth_components = self.confs['smooth_components']
        if self.confs['smooth_components'] is None:
            self._smooth_components = (set(self.confs['weights'].keys()) - set(['depth',
                                                                  'obs',
                                                                  'bound',
                                                                  'refers']))
        self.constraints = dict()
        self.constraints_val = dict()
        if self.confs['depth_constraint'] is None:
            self.constraints['depth'] = np.ones(np.prod(self.nzyx))
            self.constraints_val['depth'] = None
        else:
            self.constraints['depth'] = (self.confs['depth_constraint'].reshape(-1,1)*np.ones((1,self._nx*self._ny))).ravel()
            self.constraints_val['depth'] = 0
        if self.confs['refer_densities'] is None:
            self.constraints['refers'] = None
            self.constraints_val['refers'] = None
        else:
            self.constraints['refers'] = np.ones(self._nx*self._ny*self._nz)
            self.refer_densities = self.confs['refer_densities']
            self.n_refer = len(self.confs['refer_densities'])
        self._weights = self.confs['weights']
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
        self.optimize_log = {'parameters':[],'abic_vals':[]}

    def gen_model_name(self):
        '''generate a file name to save data of current model. The model will be
        saved in ``self.data_dir`` directory.
        '''
        self.model_name = '{}x{}x{}'.format(self._nx,self._ny,self._nz)
        self.fname = pathlib.Path(self.data_dir)/pathlib.Path(self.model_name+'.h5')

    @property
    def gpuid(self):
        ''' Which gpu card will be used. Ignored if set ``use_gpu=0``.
        '''
        return self.confs['gpuid']
    @gpuid.setter
    def gpuid(self,values):
        self.confs['gpuid'] = values

    @property
    def smooth_on(self):
        '''Which variable should be smoothed. Only ``'m'``
        is supported right now, which means smooth on density.
        '''
        return self.confs['smooth_on']
    @smooth_on.setter
    def smooth_on(self,values):
        self.confs['smooth_on'] = values

    @property
    def mode(self):
        '''How to calculate determinants. Could be 'naive' or 'walsh'. '''
        return self.confs['mode']
    @mode.setter
    def mode(self,values):
        self.confs['mode'] = values

    @property
    def data_dir(self):
        '''Specify a path to save model data.
        '''
        return self.confs['data_dir']
    @data_dir.setter
    def data_dir(self,values):
        self.confs['data_dir'] = values

    @property
    def source_volume(self):
        ''' The extent of source volume, in the form of ``[xmin,xmax,ymin,ymax,zmin,zmax]``.
        '''
        return self.confs['source_volume']
    @source_volume.setter
    def source_volume(self,value):
        self._source_volume = value
        self.confs['source_volume'] = value
        self.gen_mesh()

    @property
    def subtract_mean(self):
        return self.confs['subtract_mean']
    @subtract_mean.setter
    def subtract_mean(self,values):
        self.confs['subtract_mean'] = values

    @property
    def weights(self):
        ''' inverse variance of each distributions.
        '''
        return self._weights
    @weights.setter
    def weights(self,values):
        self._weights = values
        self.confs['weights'] = values
        if not self.kernel_op is None:
            self.kernel_op.weights = self._weights

    @property
    def optimize_keys(self):
        ''' inverse variance of each distributions.
        '''
        return self.confs['optimize_keys']
    @optimize_keys.setter
    def optimize_keys(self,values):
        self.confs['optimize_keys'] = values

    @property
    def smooth_components(self):
        ''' partial derivatives used as smooth components.
        Example: ``'dxx'`` means :math:`\frac{\partial^2 m}{\partial x^2}`
        '''
        return self._smooth_components
    @smooth_components.setter
    def smooth_components(self,values):
        self._smooth_components = values
        self.confs['smooth_components'] = self._smooth_components


    @property
    def refer_densities(self):
        '''reference density. The length of this vector should match the length
        of model density.
        '''
        tmp = []
        for density in self.constraints_val['refers']:
            tmp.append(density.reshape(self._nz,self._ny,self._nx))
        return tmp
    @refer_densities.setter
    def refer_densities(self,value):
        tmp = []
        for density in value:
            tmp.append(density.ravel())
        self.constraints_val['refers'] = tmp

    @property
    def nzyx(self):
        '''model dimension, with the form of ``[nz,ny,nx]``
        '''
        return self.confs['nzyx']
    @nzyx.setter
    def nzyx(self,values):
        self.confs['nzyx'] = values
        self._nz,self._ny,self._nx = values
        self.nobsx = values[2]
        self.nobsy = values[1]
        self.gen_model_name()
        if not self.constraints['depth'] is None:
            self.constraints['depth'] = self.constraints['depth'].reshape(self._nz,-1)[:,0]*np.ones((1,self._nx*self._ny))
            self.constraints['depth'] = self.constraints['depth'].ravel()
        self.constraints['refers'] = np.ones(self._nx*self._ny*self._nz)

    @property
    def nx(self):
        ''' Number of cells along x-axis.
        '''
        return self._nx
    @nx.setter
    def nx(self,value):
        self._nx = value
        self.nobsx = self._nx
        self.gen_model_name()
        if not self.constraints['depth'] is None:
            self.constraints['depth'] = self.constraints['depth'].reshape(self._nz,-1)[:,0]*np.ones((1,self._nx*self._ny))
            self.constraints['depth'] = self.constraints['depth'].ravel()
        self.constraints['refers'] = np.ones(self._nx*self._ny*self._nz)

    @property
    def ny(self):
        ''' Number of cells along y-axis.
        '''
        return self._ny
    @ny.setter
    def ny(self,value):
        self._ny = value
        self.nobsy = self._ny
        self.gen_model_name()
        if not self.constraints['depth'] is None:
            self.constraints['depth'] = self.constraints['depth'].reshape(self._nz,-1)[:,0]*np.ones((1,self._nx*self._ny))
            self.constraints['depth'] = self.constraints['depth'].ravel()
        self.constraints['refers'] = np.ones(self._nx*self._ny*self._nz)

    @property
    def nz(self):
        ''' Number of cells along z-axis.
        '''
        return self._nz
    @nz.setter
    def nz(self,value):
        self._nz = value
        self.gen_model_name()
        self.constraints['refers'] = np.ones(self._nx*self._ny*self._nz)
        print("Warning: nz changed. \nDon't forget setting depth constraints.")

    @property
    def model_density(self):
        ''' This vector is used for calculating gravity field, i.e. forward calculating.
        '''
        return(self._model_density.reshape(self.nz,self.ny,self.nx))
    @model_density.setter
    def model_density(self,value):
        self._model_density = value.ravel()
        self.confs['model_density'] = self._model_density

    def gen_mesh(self,height = -1):
        ''' Generate mesh of the model.
        Args:
            height (float): height of the observations.
        '''
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
        '''generate walsh matrix in the order of sequence2.
        '''
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
        walsh_matrix = walsh_matrix.astype(np.float64)
        step = self._nx*self._ny*self._nz//4
        components = ['0','1','2','3']
        with h5py.File(self.fname,mode='a') as f:
            fgroup = f.create_group('walsh_matrix')
            for i in range(4):
                fgroup.create_dataset(components[i],data=walsh_matrix[i*step:(i+1)*step,:])

    def gen_kernel(self, process = 1):
        '''generate kernel matrix. Because it is multilevel toeplitz matrix, only
        the first row are needed (i.e. all source cell contribution to the first
        observation position).
        .. note:: The coordinate system of the input parameters is to be
            x -> North, y -> East and z -> **DOWN**.
        .. note:: All input values in **SI** units(!) and output in **mGal**!
        '''
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
                                         refer_constraints=self.constraints['refers'],
                                         weights=self._weights)

    def load_kernel(self,fname):
        '''load kernel matrix from file. Only the first row are needed.
        File format follows numpy's savetxt.
        '''
        try:
            self.kernel0 = np.loadtxt(fname)
        except OSError:
            fname = pathlib.Path(self.data_dir)/pathlib.Path(fname)
            self.kernel0 = np.loadtxt(fname)
        self.kernel0 = np.array(self.kernel0).reshape(self.nz,self.ny,self.nx)
        self.kernel_op = AbicLSQOperator(self.kernel0,
                                         depth_constraint=self.constraints['depth'],
                                         smooth_components=self._smooth_components,
                                         refer_constraints=self.constraints['refers'],
                                         weights=self._weights)

    def _diagvec(self,vec=None,diag=None):
        if vec.ndim == 1:
            return vec * diag
        else:
            return  vec * diag.reshape(1,-1)

    @timeit
    def walsh_transform(self,keys=None):
        '''walsh transform of kernel matrix and constraint matrices.
        '''
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
        is_stored['refers'] = True
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
        '''Diagonal of the depth constraint matrix.
        One real number for each layer. Stored in an vector.
        '''
        return(self.constraints['depth'].reshape(self._nz,-1)[:,0])
    @depth_constraint.setter
    def depth_constraint(self,value):
        self.constraints['depth'] = (value.reshape(-1,1)*np.ones((1,self._nx*self._ny))).ravel()

    @timeit
    def forward(self,model_density=None):
        ''' Calculate gravity field from model_density.
        Args:
            model_density (np.Array): densities of each model cell. Reshaped from
            (nz,ny,nx) to (nz*ny*nx)
        '''
        if model_density is None:
            model_density = self._model_density
        else:
            model_density = model_density.ravel()
        self.obs_data = self.kernel_op.gtoep.matvec(model_density)

    def _gen_rhs(self):
        '''generate right hand side of the least square equation of min_u.
        '''
        self.rhs = self._weights['obs']*self.kernel_op.gtoep.rmatvec(self.obs_data)
        if 'refers' in self._weights.keys():
            for i,refer_weight in enumerate(self._weights['refers']):
                v = self.constraints_val['refers'][i].ravel()
                if 'depth' in self._weights.keys():
                    v = self.constraints['depth']*v
                self.rhs += (refer_weight
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
    def do_linear_solve(self,tol=1.0e-5):
        ''' solve the least square equation of min_u.
        Args:
            tol (float): tol of CG algorithm.
        '''
        self.do_linear_solve_quiet(tol=tol)

    def do_linear_solve_quiet(self,tol=1.0e-5):
        ''' solve the least square equation of min_u with minimum of message printed.
        Args:
            tol (float): tol of CG algorithm.
        '''
        self._gen_rhs()
        if self.subtract_mean:
            sum_obs = np.sum(self.obs_data)
            tmp_b = np.zeros(len(self.rhs)+1) #np.zeros(len(self.rhs+1)) #chenshi
            tmp_b[:-1] = self.rhs
            tmp_b[-1] = sum_obs
            tmp_op = AbicLSQOperator2(self.kernel_op)
            self.solution = spsparse.linalg.cg(tmp_op,tmp_b,tol=tol)[0]
        else:
            self.solution = spsparse.linalg.cg(self.kernel_op,self.rhs,tol=tol)[0]

    @timeit
    def calc_u(self,solved=False,x=None):
        '''calc min_u value.
        Args:
            solved (bool): whether or not the least square equation solved already.
            x (array): use this x to calculate value of U instead of doing optimization.
        '''
        return self.calc_u_quiet(solved,x)

    @timeit
    def calc_min_u(self,solved=False,x=None):
        '''Another name of calc_u. We keep this function for backward compatability.
        '''
        return self.calc_u_quiet(solved,x)

    def calc_u_quiet(self,solved=False,x=None):
        '''calc min_u value with minimum message printed out.
        Args:
            solved (bool): whether or not the least square equation solved already.
            x (array): use this x to calculate value of U instead of doing optimization.
        '''
        if x is None:
            if not solved:
                self.do_linear_solve_quiet()
            x = self.solution
        self.min_u_val = self._weights['obs']*np.linalg.norm(self.kernel_op.gtoep.matvec(x) - self.obs_data)**2
        if ('refers' in self._weights.keys()) and (self.smooth_on == 'm-m0'):
            v = x - self.constraints_val['refer']
        else:
            v = x
        if 'depth' in self._weights.keys():
            v = self.constraints['depth']*v
        for key in self._smooth_components:
            tmp2 = self.smop.derivation(v.reshape(self._nz,self._ny,self._nx),
                                        component=key)
            self.min_u_val += self._weights[key]*np.linalg.norm(tmp2.ravel())**2
        if 'refers' in self._weights.keys():
            for i,refer_density in enumerate(self.constraints_val['refers']):
                v = x - refer_density
                if 'depth' in self._weights.keys():
                    v = self.constraints['depth']*v
                self.min_u_val += self._weights['refers'][i] *np.linalg.norm(v)**2
        return self.min_u_val

    def jac_u(self,x=None):
        ''' jacobian of the function U.
        '''
        res = self.kernel_op.matvec(x) - self.rhs
        return 2.*res

    def hessp_u(self,x,v):
        '''Hessian of the function U.
        '''
        res = self.kernel_op.matvec(v)
        return 2.*res

    @timeit
    def bound_optimize(self,x0=None):
        '''optimize function U using boundary constraint.
        '''
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
        '''Minimize the 1-norm instead 2-norm of the model equation. This function
           is used to the target function.
        '''
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
        '''jacobian of the lasso target function.
        '''
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
        '''hessian of the lasso target function.
        '''
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
        '''optimize the lasso function.
        '''
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
        '''calculate the residual information including residuals and stds.
        The result of this function are stored in ``self.residuals`` and ``self.stds``.
        '''
        self.residuals = dict()
        self.stds = dict()
        self.residuals['obs'] = np.linalg.norm(self.kernel_op.gtoep.matvec(self.solution)-self.obs_data)**2
        self.stds['obs'] = np.std(self.kernel_op.gtoep.matvec(self.solution)-self.obs_data)
        for key in self._smooth_components:
            tmp2 = self.solution.reshape(self._nz,self._ny,self._nx)
            if ('refers' in self.constraints_val.keys()) and (self.smooth_on == 'm-m0'):
                tmp2 -= self.constraints_val['refer'].reshape(self._nz,self._ny,self._nx)
            tmp2 = self.smop.derivation(tmp2,component=key)
            self.residuals[key] = np.linalg.norm(tmp2.ravel())**2
            self.stds[key] = np.std(tmp2.ravel())
        if 'refers' in self.constraints_val.keys():
            self.residuals['refers'] = []
            self.stds['refers'] = []
            for i,refer_density in self.constraints_val['refers']:
                self.residuals['refers'].append(np.linalg.norm(self.solution.ravel()-refer_density.ravel())**2)
                self.stds['refers'].append(np.std(self.solution.ravel()-refer_density.ravel()))

    def calc_log_prior_total_det_naive(self):
        '''calculate the determinant of prior distribution and joint distribution
        with minimum message printed out.
        '''
        self.log_prior_det_val = 0
        self.log_total_det_val = 0
        prior_eigs = np.zeros(self._nx*self._ny*self._nz)
        total_eigs = np.zeros(self._nx*self._ny*self._nz)
        tmp_mat = np.zeros((self._nz*self._ny*self._nx,self._nz*self._ny*self._nx))
        for dxyz_name in self._smooth_components:
            tmp_mat += self._weights[dxyz_name]*self.matrices[dxyz_name]
        if 'depth' in self._weights.keys():
            tmp_mat = self.constraints['depth'].reshape(-1,1)*self.constraints['depth'].reshape(1,-1)*tmp_mat
        prior_eigs = np.linalg.svd(tmp_mat,compute_uv=False)
        eps = prior_eigs.max() * len(prior_eigs) * np.finfo(np.float64).eps
        eigs = prior_eigs[prior_eigs>eps]
        self.log_prior_det_val = np.sum(np.log(eigs))

        tmp_mat += self._weights['obs']*self.matrices['obs']
        if 'refers' in self._weights.keys():
            tmp_mat += sum(self._weights['refers'])*np.diag(self.constraints['depth'])**2
        uu,total_eigs,vv = np.linalg.svd(tmp_mat,compute_uv=True)
        self.log_total_det_val = np.sum(np.log(total_eigs))
        self._gen_rhs()
        self.solution = np.zeros(np.prod(self.nzyx))
        self.solution = (vv.T @ ((1./total_eigs).ravel() * (uu.T @ self.rhs)))

        self.eigs = {'prior':prior_eigs,'total':total_eigs}
        return self.log_prior_det_val,self.log_total_det_val

    def calc_log_prior_total_det_walsh(self):
        '''calculate the determinant of prior distribution and joint distribution
        with minimum message printed out.
        '''
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
            self._gen_rhs()
            self.solution = np.zeros(np.prod(self.nzyx))
            for i_b,block in enumerate(blocks):
                tmp_block = np.zeros((step,step))
                walsh_group = f['walsh_matrix']
                part_walsh = walsh_group[blocks[i_b]][:]
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

                if use_gpu > 0:
                    import cupy as cp
                    with cp.cuda.Device(self.gpu_id):
                        tmp_block_gpu = cp.asarray(tmp_block,dtype=np.float64)
                        eigs = cp.linalg.svd(tmp_block_gpu,compute_uv=False)
                        prior_eigs[i_b*step:(i_b+1)*step] = cp.asnumpy(eigs)
                        eps = eigs.max() * len(eigs) * np.finfo(np.float64).eps
                        eigs = eigs[eigs>eps]
                        self.log_prior_det_val += cp.asnumpy(cp.sum(cp.log(eigs)))
                        tmp_block_gpu = None
                        eigs = None
                        free_gpu()
                    tmp_block += self._weights['obs']*f['kernel'][block][:]
                    if 'refers' in self._weights.keys():
                        tmp_multi_small = depth_walsh.T@depth_walsh
                        for i in range(step//self._nz):
                            tmp_block[i*self._nz:(i+1)*self._nz,
                                      i*self._nz:(i+1)*self._nz] += sum(self._weights['refers'])*tmp_multi_small
                    with cp.cuda.Device(self.gpu_id):
                        tmp_block_gpu = cp.asarray(tmp_block,dtype=np.float64)
                        eigs = cp.linalg.svd(tmp_block_gpu,compute_uv=False)
                        total_eigs[i_b*step:(i_b+1)*step] = cp.asnumpy(eigs)
                        #eigs = eigs[eigs>1.0e-12]
                        self.log_total_det_val += cp.asnumpy(cp.sum(cp.log(eigs)))
                        tmp_block_gpu = None
                        eigs = None
                        free_gpu()
                else:
                    tmp_block_gpu = np.asarray(tmp_block,dtype=np.float64)
                    eigs = np.linalg.svd(tmp_block_gpu,compute_uv=False)
                    prior_eigs[i_b*step:(i_b+1)*step] = eigs
                    eps = eigs.max() * len(eigs) * np.finfo(np.float64).eps
                    eigs = eigs[eigs>eps]
                    self.log_prior_det_val += np.sum(np.log(eigs))
                    tmp_block_gpu = None
                    eigs = None
                    tmp_block += self._weights['obs']*f['kernel'][block][:]
                    if 'refers' in self._weights.keys():
                        tmp_multi_small = depth_walsh.T@depth_walsh
                        #tmp_multi_small = np.eye(tmp_multi_small.shape[0])
                        for i in range(step//self._nz):
                            tmp_block[i*self._nz:(i+1)*self._nz,
                                      i*self._nz:(i+1)*self._nz] += sum(self._weights['refers'])*tmp_multi_small
                    tmp_block_gpu = np.asarray(tmp_block,dtype=np.float64)
                    uu,eigs,vv = np.linalg.svd(tmp_block_gpu,compute_uv=True)
                    total_eigs[i_b*step:(i_b+1)*step] = eigs
                    #eigs = eigs[eigs>1.0e-12]
                    self.log_total_det_val += np.sum(np.log(eigs))
                    tmp_block_gpu = None
                    self.solution += part_walsh.T @ (vv.T @ ((1./eigs).ravel() * (uu.T @ (part_walsh @ self.rhs))))
                    eigs = None


        if use_gpu > 0:
            self.log_prior_det_val = cp.asnumpy(self.log_prior_det_val)
            self.log_total_det_val = cp.asnumpy(self.log_total_det_val)
        else:
            self.log_prior_det_val = self.log_prior_det_val
            self.log_total_det_val = self.log_total_det_val

        self.eigs = {'prior':prior_eigs,'total':total_eigs}
        return self.log_prior_det_val,self.log_total_det_val

    def calc_log_prior_total_det_quiet(self,mode=None):
        '''calculate the determinant of prior distribution and joint distribution
        with minimum message printed out.
        '''
        if mode is None:
            mode = self.mode
        if mode == 'walsh':
            self.calc_log_prior_total_det_walsh()
        elif mode == 'naive':
            self.calc_log_prior_total_det_naive()
        else:
            raise ValueError('mode={} is not implemented!!'.format(mode))

    def prepare_det(self,mode=None):
        if mode is None:
            mode = self.mode
        if mode == 'walsh':
            self.walsh_transform()
        elif mode == 'naive':
            self.matrices = dict()
            n_cell = self._nx*self._ny*self._nz
            self.matrices['G']  = self.kernel_op.gtoep.matvec(np.eye(self._nx*self._ny*self._nz)).T
            self.matrices['obs'] = self.kernel_op.gtoep.rmatvec(self.kernel_op.gtoep.matvec(np.eye(n_cell)))
            for key in self._smooth_components:
                self.matrices[key] = self.smop.rderivation(self.smop.derivation(np.eye(n_cell).reshape(-1,self._nz,self._ny,self._nx),key),key).reshape(n_cell,n_cell)
        else:
            raise ValueError('mode={} is not implemented!!'.format(mode))

    @timeit
    def calc_log_prior_total_det(self,mode=None):
        '''calculate the determinant of prior distribution and joint distribution.
        '''
        return self.calc_log_prior_total_det_quiet(mode=mode)

    def calc_log_obs_det_quiet(self):
        '''calculate the determinant of observation's distribution with minimum
        message printed out.
        '''
        self.log_obs_det_val = (4*sum(np.log(self.constraints['depth']))
                                +np.prod(self.nzyx)*np.sum(np.log(self.weights['refers']))
                                 +np.log(self._weights['obs'])*len(self.obs_data))
        return self.log_obs_det_val

    @timeit
    def calc_log_obs_det(self):
        '''calculate the determinant of observation's distribution message printed out.
        '''
        return self.calc_log_obs_det_quiet()

    @timeit
    def calc_abic(self,mode=None):
        '''calculate abic value: -log_prior_det_value+log_total_det-log_obs_det+min_u'''
        self.calc_log_obs_det()
        self.calc_log_prior_total_det(mode=mode)
        self.calc_u(solved=True)
        self.abic_val = (self.log_total_det_val
                         + self.min_u_val
                         - self.log_prior_det_val
                         - self.log_obs_det_val)
        return self.abic_val

    def calc_abic_quiet(self,mode=None):
        '''calculate abic value: -log_prior_det_value+log_total_det-log_obs_det+min_u'''
        self.calc_log_obs_det_quiet()
        self.calc_log_prior_total_det_quiet(mode=mode)
        self.calc_u_quiet(solved=True)
        self.abic_val = (self.log_total_det_val
                         + self.min_u_val
                         - self.log_prior_det_val
                         - self.log_obs_det_val)
        self.optimize_log['parameters'].append(self.weights)
        self.optimize_log['abic_vals'].append(self.abic_val)
        return self.abic_val

    def _abic_optimize_exp(self,mode=None,**opt_args):
        '''optimize the abic value. Use the log and exp trick to constrain weights
        to be positive.
        '''
        #optimize_keys = list(set(self._weights.keys())-set(['depth']))
        if self.confs['optimize_keys'] is None:
            optimize_keys = list(self._weights.keys())
        else:
            optimize_keys = self.confs['optimize_keys']
        optimize_keys_1 = list(set(optimize_keys)-set(['refers']))
        def abic_target(x):
            tmp_weights = dict()
            for i,key in enumerate(optimize_keys_1):
                tmp_weights[key] = np.exp(x[i])
            if 'refers' in optimize_keys:
                tmp_weights['refers'] = list(np.exp(x[-len(self.weights['refers']):]))
            self.weights = {**self.weights,**tmp_weights}
            self.calc_abic_quiet(mode=mode)
            return self.abic_val
        n_weights = len(optimize_keys)
        if 'refers' in optimize_keys:
            n_weights += len(self.weights['refers']) - 1
        x0 = np.zeros(n_weights)
        for i,key in enumerate(optimize_keys_1):
            x0[i] = np.log(self._weights[key])
        if 'refers' in optimize_keys:
            for i,refer_weight in enumerate(self._weights['refers']):
                x0[len(optimize_keys_1)+i] = np.log(refer_weight)
        self.abic_optimize_summary = minimize(abic_target,
                                              x0,
                                              **opt_args)

    def _abic_optimize_bound(self):
        '''optimize the abic value. Use scipy's COBYLA algorithm to constrain weights
        to be positive.
        '''
        if self.confs['optimize_keys'] is None:
            optimize_keys = list(self._weights.keys())
        else:
            optimize_keys = self.confs['optimize_keys']
        def abic_target(x):
            tmp_weights = dict()
            for i,key in enumerate(optimize_keys):
                tmp_weights[key] = np.exp(x[i])
            self.weights = {**self.weights,**tmp_weights}
            return self.calc_abic_quiet()
        x0 = np.zeros(len(optimize_keys))
        for i,key in enumerate(optimize_keys):
            x0[i] = self._weights[key]
        weight_constraint = LinearConstraint(np.eye(len(optimize_keys)),0.,np.inf)
        self.abic_optimize_summary = minimize(abic_target,
                                              x0,
                                              method='COBYLA',
                                              constraints=weight_constraint)

    @timeit
    def abic_optimize(self,mode=None,**opt_args):
        '''optimize the abic value. This is the interface for users.
        '''
        self._abic_optimize_exp(mode=mode,**opt_args)

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
