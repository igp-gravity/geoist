import pathlib
from datetime import datetime
from functools import wraps
import numpy as np
from scipy import linalg as splin
from scipy import sparse as spsparse
from scipy.optimize import minimize
import pandas as pd

from geoist import gridder
from geoist.pfm import prism,tesseroid
from geoist.inversion.mesh import PrismMesh,TesseroidMesh
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid

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

class SmoothOperator:
    def __init__(self):
        self.axis = {'tzyx':{'x':-1,'y':-2,'z':-3,'t':-4},
                     'txyz':{'t':-4,'x':-3,'y':-2,'z':-1}}

    def derivation(self,v,component='dx',array_order='tzyx'):
        for axis_i in component[1:]:
            slices = [slice(None)]*v.ndim
            slices[self.axis[array_order][axis_i]] = slice(-1,None,-1)
            v = np.diff(v[tuple(slices)],axis=self.axis[array_order][axis_i])[tuple(slices)]
        return v

    def rderivation(self,v,component='dx',array_order='tzyx'):
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

class AbicLSQOperatorMS:
    '''An operator doing matrix vector multiplication. The matrix is:
        $\alpha_g G^TG + \sum \alpha_i B_i^TB_i$. Where $\alpha$'s are
        weights, $G$ is kernel matrix, $B_i$'s are smooth matrix.
    Args:
        kernel_matrix (list of ndarray): kernel matrix of each survey
        ns (int): number of surveys
        nzyx (tuple of int): number of cells along z-axis, y-axis and x-axis
        smmoth_components (list of str): which components should be smoothed, acceptable
            string could be 'dx','dxx','dy','dyy','dxy','dt',...
        weights (dict): name and values of weights.
    '''
    def __init__(self,
                 kernel_matrices,
                 ns,
                 nzyx,
                 smooth_components=set(),
                 weights=None):
        self.kernel_matrices = kernel_matrices
        self.ns = ns
        self.nz,self.ny,self.nx = nzyx
        self.shape = (self.ns*self.ny*self.nx,self.ns*self.ny*self.nx)
        self.weights = weights
        self.smooth_components = smooth_components
        self.smop = SmoothOperator()
        if self.weights is None:
            self.weights = {'bound':1,'obs':1,'dx':1,'dy':1,'dt':1, 'refer':1}

    def matvec(self,v):
        tmp = np.zeros_like(v)
        v = v.reshape(self.ns,-1)
        for i,k in enumerate(self.kernel_matrices):
            tmp[i*self.ny*self.nx:(i+1)*self.ny*self.nx] = self.weights['obs'] * (k.T @ (k @ v[i]))
        if 'refer' in self.weights.keys():
            tmp += self.weights['refer']*v.ravel()
        for key in self.smooth_components:
            tmp2 = v.reshape(self.ns,self.nz,self.ny,self.nx)
            tmp2 = self.smop.derivation(tmp2,component=key)
            tmp2 = self.smop.rderivation(tmp2,component=key)
            tmp += self.weights[key]*tmp2.ravel()
        return tmp

class InvModelTS:
    '''Multiple survey gravity inversion use Abic
    Args:
        nyx (tuple of int): number of cells along y-axis and x-axis
        smooth_components (list of str): which components should be smoothed, acceptable
            string could be 'dx','dxx','dy','dyy','dxy','dt',...
        weights (dict): name and values of weights.
        optimize_weights (list of str): specify which weight should be optimized by Abic.
        source_volume (tuple of float): dimension of the underlying gravity source.
        margin (tuple of float): margin size around the source volume.
        cell_type (str): either 'prism' or 'Tesseroid'
        data_dir (str): a folder where observation data resides
    Attributes:
        ns (int): number of surveys
        nx,ny,nz (int): number of cells along x-axis, y-axis and z-axis
        source_volume (tuple of float): dimension of the underlying gravity source
        margin (tuple of float): margin size around the source volume.
        cell_type (str): either 'prism' or 'Tesseroid'
        smooth_components (list of str): which components should be smoothed
        weights (dict): name and values of weights.
        optimize_weights (list of str): specify which weight should be optimized by Abic.
        data_dir (str): a folder where observation data resides
        abic_val (float): abic value
        log_total_det_val (float): log determinate value of the total matrix
        log_obs_det_val (float): log determinate value of the observation matrix
        min_u_val (float): minimum of U
        min_density (float): minimum density
        max_density (float): maximum density
        kernel_matrices (list of ndarray): kernels of each survey.
        rhs (ndarray): right hand side when solve min_u
        solution (ndarray): model densityies solve min_u
        abic_log (dict): save the weights and abic value during calculation
        orig_data (DataFrame): original data loaded from files.
        mesh (Mesh): mesh of the underlying density model.
        R (float): radius of the Earth (m)
    '''
    def __init__(self,
                 nyx=(20,20),
                 smooth_components=None,
                 weights=None,
                 optimize_weights=None,
                 source_volume=None,
                 margin=(0,0,0,0),
                 cell_type='prism',
                 data_dir='./data'):
        self.ny,self.nx = nyx
        self.nz = 1
        self.data_dir = pathlib.Path(data_dir)
        self.source_volume = source_volume
        self.margin= margin
        self.kernel_op = None
        self.smooth_components = smooth_components
        self._weights = weights
        self.optimize_weights = optimize_weights
        self.cell_type = cell_type.lower()
        self.smop = SmoothOperator()
        self.kernel_matrices = []
        self.abic_val = 0
        self.log_total_det_val = 0
        self.log_prior_det_val = 0
        self.log_obs_det_val = 0
        self.min_u_val = 0
        self.min_density = -1.0e4
        self.max_density = 1.0e4
        self._btb_exist = False
        self._gtg_exist = False
        self.abic_log = {'weights':[],'abic_val':[]}
        self._abic_iter = 0
        self.R = 6371000
        self._keymap = None
        self.orig_data = None

    def load_data(self,pattern='*.txt',names=['lon','lat','g'],**kwargs):
        ''' load all data inside data dir, result saved in self.orig_data.
        Args:
            pattern (str): filename pattern.
            names (list of str): name of the observation must be 'g', coordinates
                must be 'lon','lat' or 'x','y'
            kwargs : same as pd.read_csv
        '''
        data_files = sorted(self.data_dir.glob(pattern))
        self.ns = len(data_files)
        if not kwargs:
            kwargs = dict()
            kwargs['delim_whitespace'] = True
        if self.cell_type == 'prism':
            self._keymap = ['x','y']
        elif self.cell_type == 'tesseroid':
            self._keymap = ['lon','lat']
        orig_data = []
        for i,data_file in enumerate(data_files):
            df = pd.read_csv(data_file,names=names,**kwargs)
            i_survey = np.zeros(len(df),dtype=int) + i
            df['i_survey'] = i_survey
            orig_data.append(df)
        self.orig_data = pd.concat(orig_data)
        self.orig_data['z'] = 0.0

    def set_refer(self,refer_density):
        self.refer_density = refer_density

    def deg2xy(self):
        dlon = self.orig_data.groupby('i_survey')['lon'].apply(lambda x: x-x.mean())
        dlat = self.orig_data.groupby('i_survey')['lat'].apply(lambda x: x-x.mean())
        x = dlat*self.R*np.pi/180.
        y = dlon*self.R*np.cos(self.orig_data['lat']*np.pi/180.)*np.pi/180.
        self.orig_data['x'] = x
        self.orig_data['y'] = y

    @property
    def weights(self):
        return self._weights
    @weights.setter
    def weights(self,values):
        self._weights = values
        if not self.kernel_op is None:
            self.kernel_op.weights = self._weights
    @property
    def nx(self):
        return self._nx
    @nx.setter
    def nx(self,value):
        self._nx = value
        self._btb_exist = False
        self._gtg_exist = False

    @property
    def ny(self):
        return self._ny
    @ny.setter
    def ny(self,value):
        self._ny = value
        self._btb_exist = False
        self._gtg_exist = False

    @property
    def smooth_components(self):
        return self._smooth_components
    @smooth_components.setter
    def smooth_components(self,values):
        self._smooth_components = values
        if not self.kernel_op is None:
            self.kernel_op.smooth_components = self._smooth_components

    def _gen_source_volume(self,source_volume=None,margin=None):
        if source_volume is None:
            source_volume = list(self.source_volume)
        fun_list = [min,max,min,max]
        key_ind = [self._keymap[0],self._keymap[0],self._keymap[1],self._keymap[1]]
        for i in range(4):
            if source_volume[i] is None:
                source_volume[i] = fun_list[i](self.orig_data[key_ind[i]])
                source_volume[i] += margin[i]*(-1)**(i+1)
        self.source_volume = source_volume

    def gen_mesh(self,source_volume=None,margin=None):
        shape = (self.nz, self.ny, self.nx)
        if margin is None:
            margin = self.margin
        self._gen_source_volume(source_volume,margin)
        if self.cell_type =='prism':
            self.mesh = PrismMesh(self.source_volume, shape)
        elif self.cell_type =='tesseroid':
            self.mesh = TesseroidMesh(self.source_volume, shape)
        else:
            raise ValueError('cell_type must be \'prism\' or \'tesseroid\'!!')
        density = np.ones(shape)*1.0e3
        self.mesh.addprop('density', density.ravel())

    def gen_obs_grid(self,height=-1):
        """ generate obs grid
         """
        obs_area = (self.source_volume[0]+0.5*self.mesh.dims[0],
                         self.source_volume[1]-0.5*self.mesh.dims[0],
                         self.source_volume[2]+0.5*self.mesh.dims[1],
                         self.source_volume[3]-0.5*self.mesh.dims[1])
        obs_shape = (self.nx, self.ny)
        return gridder.regular(obs_area, obs_shape, z=height)

    def gen_kernel(self):
        self.kernel_matrices = []
        if self.cell_type == 'prism':
            for i in range(self.ns):
                xp = self.orig_data[self.orig_data['i_survey']==i]['x'].values
                yp = self.orig_data[self.orig_data['i_survey']==i]['y'].values
                zp = self.orig_data[self.orig_data['i_survey']==i]['z'].values
                kernel0 = np.zeros((len(xp),self.ny*self.nx))
                for j,cell in enumerate(self.mesh):
                    kernel0[:,j] = prism.gz(xp,yp,zp,[cell])
                self.kernel_matrices.append(np.array(kernel0))
        elif self.cell_type == 'tesseroid':
            for i in range(self.ns):
                xp = self.orig_data[self.orig_data['i_survey']==i]['lon'].values
                yp = self.orig_data[self.orig_data['i_survey']==i]['lat'].values
                zp = self.orig_data[self.orig_data['i_survey']==i]['z'].values
                kernel0 = np.zeros((len(xp),self.ny*self.nx))
                for j,cell in enumerate(self.mesh):
                    kernel0[:,j] = tesseroid.gz(xp,yp,zp,[cell])
                self.kernel_matrices.append(np.array(kernel0))
        else:
            raise ValueError('cell_type must be \'prism\' or \'tesseroid\'!!')
        self.kernel_op = AbicLSQOperatorMS(self.kernel_matrices,
                                           ns=self.ns,
                                           nzyx=(self.nz,self.ny,self.nx),
                                         smooth_components=self._smooth_components,
                                         weights=self._weights)
        self._gtg_exist = False

    def _diagvec(self,vec=None,diag=None):
        if vec.ndim == 1:
            return vec * diag
        else:
            return  vec * diag.reshape(1,-1)

    @timeit
    def forward(self,model_density=None):
        obs_g = []
        if model_density is None:
            x = self.solution.reshape(self.ns,-1)
        else:
            x = model_density.reshape(self.ns,-1)
        for i,k in enumerate(self.kernel_matrices):
            obs_g.append((k @ x[i]).ravel())
        return np.hstack(obs_g)

    def _gen_rhs(self):
        self.rhs = np.zeros(self.ns*self.ny*self.nx)
        for i,k in enumerate(self.kernel_matrices):
            s = i*self.ny*self.nx
            e = (i+1)*self.ny*self.nx
            g = self.orig_data[self.orig_data['i_survey']==i]['g'].values
            self.rhs[s:e] = self._weights['obs'] * (k.T @ g)
        if 'refer' in self._weights.keys():
            self.rhs += (self._weights['refer']*self.refer_density.ravel())

    def _gen_btb(self):
        if self._btb_exist:
            return
        self._btb = dict()
        for key in self.smooth_components:
            if 't' in key:
                tmp = np.eye(self.ns*self.nx*self.ny).reshape(-1,self.ns,self.nz,self.ny,self.nx)
                nrow = self.ns*self.ny*self.nx
            else:
                tmp = np.eye(self.nx*self.ny).reshape(-1,self.ny,self.nx)
                nrow = self.nx*self.ny
            self._btb[key] = self.smop.rderivation(self.smop.derivation(tmp,component=key),
                                                   component=key).reshape(nrow,-1)
        self._btb_exist = True

    def _gen_gtg(self):
        if self._gtg_exist:
            return
        self._gtg = []
        for k in self.kernel_matrices:
            self._gtg.append(k.T @ k)
        self._gtg_exist = True

    @timeit
    def do_linear_solve(self):
        self.do_linear_solve_quiet()

    def do_linear_solve_quiet(self):
        self._gen_rhs()
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
        v = x.reshape(self.ns,-1)
        self.min_u_val = 0.
        for i,k in enumerate(self.kernel_matrices):
            g = self.orig_data[self.orig_data['i_survey']==i]['g'].values
            self.min_u_val += self._weights['obs']*np.linalg.norm(k @ v[i] - g)**2
        for key in self._smooth_components:
            tmp2 = self.smop.derivation(v.reshape(self.ns,self.nz,self.ny,self.nx),
                                        component=key)
            self.min_u_val += self._weights[key]*np.linalg.norm(tmp2.ravel())**2
        if 'refer' in self._weights.keys():
            v = x - self.refer_density.ravel()
            self.min_u_val += self._weights['refer'] * np.linalg.norm(v)**2
        return self.min_u_val

    def calc_res(self):
        self.residuals = dict()
        self.stds = dict()
        x = self.solution.reshape(self.ns,-1)
        res = []
        for i,k in enumerate(self.kernel_matrices):
            g = self.orig_data[self.orig_data['i_survey']==i]['g'].values
            res.append((k@x[i] - g).ravel())
        res = np.hstack(res)
        self.residuals['obs'] = np.linalg.norm(res)**2
        self.stds['obs'] = np.std(res)
        for key in self._smooth_components:
            tmp2 = self.solution.reshape(self.ns,self.nz,self._ny,self._nx)
            tmp2 = self.smop.derivation(tmp2,component=key)
            self.residuals[key] = np.linalg.norm(tmp2.ravel())**2
            self.stds[key] = np.std(tmp2.ravel())
        if 'refer' in self._weights.keys():
            self.residuals['refer'] = []
            self.stds['refer'] = []
            self.residuals['refer'].append(np.linalg.norm(self.solution.ravel()-self.refer_density.ravel())**2)
            self.stds['refer'].append(np.std(self.solution.ravel()-self.refer_density.ravel()))

    def calc_log_prior_total_det_quiet(self,precision=1.0e-6):
        self._gen_gtg()
        self._gen_btb()
        self.log_prior_det_val = 0
        self.log_total_det_val = 0
        prior_eigs = np.zeros(self._nx*self._ny*self.nz)
        total_eigs = np.zeros(self._nx*self._ny*self.nz)
        tmp_mat = np.zeros((self.ns*self.nx*self.ny,self.ns*self.nx*self.ny))
        for key in self.smooth_components:
            if 't' in key:
                tmp_mat += self.weights[key] * self._btb[key]
            else:
                for i in range(self.ns):
                    tmp_mat[i*self.nx*self.ny:(i+1)*self.nx*self.ny,
                            i*self.nx*self.ny:(i+1)*self.nx*self.ny] += self.weights[key] * self._btb[key]
        prior_eigs = np.linalg.eigvalsh(tmp_mat)
        self.log_prior_det_val = sum(np.log(prior_eigs[prior_eigs > precision]))
        for i,k in enumerate(self.kernel_matrices):
            tmp_mat[i*self.nx*self.ny:(i+1)*self.nx*self.ny,
                    i*self.nx*self.ny:(i+1)*self.nx*self.ny] += self.weights['obs'] * self._gtg[i]
        if 'refer' in self._weights.keys():
            tmp_mat += self._weights['refer']*np.eye(self.nx*self.ny*self.ns)
        total_eigs = np.linalg.eigvalsh(tmp_mat)
        self.log_total_det_val = sum(np.log(total_eigs[total_eigs > precision]))
        self.eigs = {'prior':prior_eigs,'total':total_eigs}

    @timeit
    def calc_log_prior_total_det(self,precision=1.0e-6):
        return self.calc_log_prior_total_det_quiet(precision)

    def calc_log_obs_det_quiet(self):
        if 'refer' in self._weights.keys():
            self.log_obs_det_val = (self.nx*self.ny*self.ns*np.log(self.weights['refer'])
                                    +np.log(self._weights['obs'])*len(self.orig_data))
        else:
            self.log_obs_det_val = np.log(self._weights['obs'])*len(self.orig_data)
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
        self.abic_log['weights'].append(self._weights.copy())
        self.abic_log['abic_val'].append(self.abic_val)
        return self.abic_val

    def calc_abic_quiet(self,precision=1.0e-6):
        '''-log_prior_det_value+log_total_det-log_obs_det+min_u'''
        self.calc_log_prior_total_det_quiet(precision)
        self.calc_u_quiet()
        self.calc_log_obs_det_quiet()
        self.abic_val = (self.log_total_det_val
                         + self.min_u_val
                         - self.log_prior_det_val
                         - self.log_obs_det_val)
        self.abic_log['weights'].append(self._weights.copy())
        self.abic_log['abic_val'].append(self.abic_val)
        if self._abic_iter % 10 == 0:
            print('abic value is:{}'.format(self.abic_val))
        self._abic_iter += 1
        return self.abic_val

    def _abic_optimize_exp(self,precision=1.0e-6):
        def abic_target(x):
            for i,key in enumerate(self.optimize_weights):
                self._weights[key] = np.exp(x[i])
            return self.calc_abic_quiet(precision=1.0e-6)
        x0 = np.zeros(len(self.optimize_weights))
        for i,key in enumerate(self.optimize_weights):
            x0[i] = np.log(self._weights[key])
        self.abic_optimize_summary = minimize(abic_target,
                                              x0,
                                              method='Nelder-Mead')

    @timeit
    def abic_optimize(self):
        self._abic_iter = 0
        self._abic_optimize_exp(precision=1.0e-6)

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

    def plot_density(self,density=None,surveys=None,fname=None):
        if surveys is None:
            surveys = range(self.ns)
        fig = plt.figure(figsize=(10, 10))
        nrows = int(np.ceil(np.sqrt(len(surveys))))
        grid = ImageGrid(fig, 111,
                        nrows_ncols=(nrows, nrows),
                        axes_pad=0.05,
                        cbar_mode='single',
                        cbar_location='right',
                        cbar_pad=0.1
                        )
        if density is None:
            x = self.solution.reshape(self.ns,self.ny,self.nx)
        else:
            x = density.reshape(self.ns,self.ny,self.nx)
        if self.cell_type == 'prism':
            #rint(x.shape)
            x = np.transpose(x,axes=[0,2,1]) #axis //chenshi
        for ind,i_survey in enumerate(surveys):
            grid[ind].set_axis_off()
            im = grid[ind].imshow(x[i_survey],origin='lower')
        cbar = grid.cbar_axes[0].colorbar(im)
        if fname is None:
            plt.show()
        else:
            plt.savefig(fname,dpi=150)

    def plot_field(self,field=None,surveys=None,fname=None,plot_station=True):
        if surveys is None:
            surveys = range(self.ns)
        if field is None:
            obs_g = self.orig_data['g']
        else:
            obs_g = pd.Series(field,index=self.orig_data.index)
        fig = plt.figure(figsize=(10, 10))
        if self.cell_type == 'prism':
            axis_order = ['y','x']
        elif self.cell_type == 'tesseroid':
            axis_order = ['lon','lat']
        nrows = int(np.ceil(np.sqrt(len(surveys))))
        grid = ImageGrid(fig, 111,
                        nrows_ncols=(nrows, nrows),
                        axes_pad=0.05,
                        cbar_mode='single',
                        cbar_location='right',
                        cbar_pad=0.1
                        )
        for ind,i_survey in enumerate(surveys):
            grid[ind].set_axis_off()
            tmp = self.orig_data[self.orig_data['i_survey']==ind]
            x = tmp[axis_order[0]].values
            y = tmp[axis_order[1]].values
            g = obs_g[self.orig_data['i_survey']==ind].values
            im = grid[ind].tricontourf(x, y, g, 20)
            if plot_station:
                im2 = grid[ind].scatter(x,y)
        cbar = grid.cbar_axes[0].colorbar(im)
        if fname is None:
            plt.show()
        else:
            plt.savefig(fname,dpi=150)
