# Copyright 2019 Pascal Audet
#
# This file is part of PlateFlex.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
"""

:mod:`~plateflex` defines the following ``Grid`` classes:

- :class:`~plateflex.classes.Grid`
- :class:`~plateflex.classes.TopoGrid`
- :class:`~plateflex.classes.GravGrid`
- :class:`~plateflex.classes.BougGrid`
- :class:`~plateflex.classes.FairGrid`
- :class:`~plateflex.classes.RhocGrid`
- :class:`~plateflex.classes.ZcGrid`

These classes can be initialized with a grid of topography/bathymetry or gravity
anomaly (Bouguer/Free-air) data, and contain methods for the following functionality:

- Extracting contours at some level of the grid
- Performing a wavelet transform using a Morlet wavelet
- Obtaining the wavelet scalogram from the wavelet transform
- Plotting the input grids, wavelet transform components, and scalograms

This module further contains the class :class:`~plateflex.classes.Project`, 
which itself is a container of :class:`~plateflex.classes.Grid` objects 
(at least one each of :class:`~plateflex.classes.TopoGrid` and 
:class:`~plateflex.classes.GravGrid`). Methods are available to:

- Add :class:`~plateflex.classes.Grid` or :class:`~plateflex.classes.Project` objects to the project
- Iterate over :class:`~plateflex.classes.Grid` objects
- Initialize the project
- Perform the wavelet admittance and coherence between topography (:class:`~plateflex.classes.TopoGrid` object) and gravity anomalies (:class:`~plateflex.classes.GravGrid` object)
- Plot the wavelet admnittance and coherence spectra
- Estimate model parameters at single grid cell
- Estimate model parameters at every (or decimated) grid cell
- Plot the statistics of the estimated parameters at single grid cell
- Plot the observed and predicted admittance and coherence functions at single grid cell
- Plot the final grids of model parameters

"""

# -*- coding: utf-8 -*-
import numpy as np
import time
import geoist.flex as gflex
from geoist.flex.cpwt import cpwt
from geoist.flex.cpwt import conf_cpwt as cf_w
from geoist.flex.flex import conf_flex as cf_f
from geoist.flex import plotting
from geoist.flex import estimate
import seaborn as sns
sns.set()


class Grid(object):
    """
    An object containing a 2D array of data and Cartesian coordinates specifying the
    bounding box of the array. Contains methods to calculate the wavelet transform, 
    wavelet scalogram and to plot those quantities at a specified wavenumber index.

    :type grid: :class:`~numpy.ndarray`
    :param grid: 2D array of of topography/gravity data
    :type dx: float
    :param dx: Grid spacing in the x-direction (km)
    :type dy: float
    :param dy: Grid spacing in the y-direction (km)

    Grid must be projected in km.

    .. rubric:: Default Attributes

    ``data`` : :class:`~numpy.ndarray`
        2D array of topography/gravity data (shape (`nx,ny`))
    ``dx`` : float 
        Grid spacing in the x-direction in km
    ``dy`` : float 
        Grid spacing in the y-direction in km
    ``nx`` : int 
        Number of grid cells in the x-direction
    ``ny`` : int 
        Number of grid cells in the y-direction
    ``units`` : str 
        Units of data set
    ``sg_units`` : str
        Units of power-spectral density of data set
    ``logsg_units`` : str
        Units of power-spectral density of data set in log
    ``title`` : str
        Descriptor for data set - used in title of plots
    ``ns`` int 
        Number of wavenumber samples
    ``k`` : np.ndarray 
        1D array of wavenumbers

    .. note:

        In all instances `x` indicates eastings in metres and `y` indicates northings.
        Using a grid of longitude / latitudinal coordinates (degrees) will result
        in incorrect calculations.

    .. rubric: Example

    >>> import numpy as np
    >>> from plateflex import Grid
    >>> # Create zero-valued square grid
    >>> nn = 200; dd = 10.
    >>> xmin = ymin = 0.
    >>> xmax = ymax = (nn-1)*dd
    >>> data = np.zeros((nn, nn))
    >>> grid = Grid(data, dx, dy)
    >>> grid
    <plateflex.grids.Grid object at 0x10613fe10>

    """

    def __init__(self, grid, dx, dy):

        nx, ny = grid.shape
        self.nx, self.ny = nx, ny
        self.dx, self.dy = dx, dy
        self.units = None
        self.sg_units = None
        self.logsg_units = None
        self.title = None
        self.ns, self.k = _lam2k(nx, ny, dx, dy)

        if np.any(np.isnan(np.array(grid))):
            
            print('grid contains NaN values. Performing interpolation...')
            
            from scipy.interpolate import griddata

            # Now interpolate grid where NaNs
            good = np.where(np.isfinite(grid))
            xx, yy = np.mgrid[0:nx, 0:ny]
            points = np.array([xx[good], yy[good]]).T
            grid_int = griddata(points, grid[good], (xx, yy), method='nearest')
            self.data = np.array(grid_int)

        else:
            self.data = np.array(grid)

    def make_contours(self, level=0.):
        """
        This function returns the contours as a List of coordinate positions
        at one given level - run this more than once with different levels 
        if desired

        :type level: float
        :param level: Level (z value) of grid to contour

        :return:
            contours: List of contours with coordinate positions

        """
        try:
            from skimage import measure
        except:
            raise(Exception("Package 'scikit-image' not available to make contours."))

        return measure.find_contours(self.data, level)

    def plot(self, mask=None, title=None, save=None, clabel=None, contours=None, **kwargs):

        if title is not None:
            title = title
        else:
            title = self.title
        if clabel is not None:
            clabel = clabel
        else:
            clabel = self.units

        plotting.plot_real_grid(self.data, title=title, mask=mask, save=save, \
            clabel=clabel, contours=contours, **kwargs)

    def wlet_transform(self):
        """
        This method uses the module :mod:`~plateflex.cpwt.cpwt` to calculate 
        the wavelet transform of the grid. By default the method mirrors the grid edges,
        tapers them (10 points on each edge) and pads the array with zeros to the next 
        power of 2. The wavelet transform is stored as an attribute of the object.

        .. rubric:: Additional Attributes

        ``wl_trans`` : :class:`~numpy.ndarray`
            Wavelet transform of the grid (shape (`nx,ny,na,ns`))
        
        .. rubric:: Example

        >>> import numpy as np
        >>> from plateflex import Grid
        >>> # Create zero-valued square array
        >>> nn = 200; dd = 10.
        >>> xmin = ymin = 0.
        >>> xmax = ymax = (nn-1)*dd
        >>> data = np.zeros((nn, nn))
        >>> grid = Grid(data, dx, dy)
        >>> grid.wlet_transform()
         #loops = 13:  1  2  3  4  5  6  7  8  9 10 11 12 13
        >>> grid.wl_trans.shape
        (100, 100, 11, 13)

        """

        # Calculate wavelet transform
        wl_trans = cpwt.wlet_transform(self.data, self.dx, self.dy, self.k)

        # Save coefficients as attribute
        self.wl_trans = wl_trans

        return

    def wlet_scalogram(self):
        """
        This method uses the module :mod:`~plateflex.cpwt.cpwt` to calculate 
        the wavelet scalogram of the grid. 

        .. note::

            If no ``wl_trans`` attribute is found, the method automatically calculates
            the wavelet transform first.

        The wavelet scalogram is stored as an attribute of the object.

        .. rubric:: Additional Attributes

        ``wl_sg`` : :class:`~numpy.ndarray`
            Wavelet scalogram of the grid (shape (`nx,ny,na,ns`))

        .. rubric:: Examples

        >>> import numpy as np
        >>> from plateflex import Grid
 
        >>> # Create random-valued square grid
        >>> nn = 200; dd = 10.
        >>> xmin = ymin = 0.
        >>> xmax = ymax = (nn-1)*dd
        >>> # set random seed
        >>> np.random.seed(0)
        >>> data = np.random.randn(nn, nn)
        >>> grid1 = Grid(data, dx, dy)
        >>> grid1.wlet_scalogram()
         #loops = 17:  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17
        >>> grid1.wl_sg.shape
        (200, 200, 17)

        >>> # Perform wavelet transform first
        >>> grid2 = Grid(data, dx, dy)
        >>> grid2.wlet_transform()
         #loops = 17:  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17
        >>> grid2.wlet_scalogram()
        >>> np.allclose(grid1.wl_sg, grid2.wl_sg)
        True
       """

        try:
            wl_sg, ewl_sg = cpwt.wlet_scalogram(self.wl_trans)
        except:
            self.wlet_transform()
            wl_sg, ewl_sg = cpwt.wlet_scalogram(self.wl_trans)

        # Save as attributes
        self.wl_sg = wl_sg
        self.ewl_sg = ewl_sg

        return

    def plot_transform(self, kindex=None, aindex=None, log=False, mask=None, \
        save=None, clabel=None, contours=None, **kwargs):
        """
        This method plots the real and imaginary components of the wavelet transform of a
        :class:`~plateflex.classes.Grid` object at wavenumber and angle indices (int). 
        Raises ``Exception`` for the cases where:

        - no wavenumber OR angle index is specified (kindex and aindex)
        - wavenumber index is lower than 0 or larger than self.ns
        - angle index is lower than 0 or larger than 11 (hardcoded)

        .. note::

            If no ``wl_trans`` attribute is found, the method automatically calculates
            the wavelet transform first.

        """

        if kindex is None or aindex is None:
            raise(Exception('Specify index of wavenumber and angle to plot the transform'))

        if kindex>self.ns-1 or kindex<0:
            raise(Exception('Invalid index: should be between 0 and '+str(self.ns)))

        if aindex>10 or aindex<0:
            raise(Exception('Invalid index: should be between 0 and 10'))

        try:
            rdata = np.real(self.wl_trans[:,:,aindex,kindex])
            idata = np.imag(self.wl_trans[:,:,aindex,kindex])
        except:
            print('Calculating the transform first')
            self.wlet_transform()
            rdata = np.real(self.wl_trans[:,:,aindex,kindex])
            idata = np.imag(self.wl_trans[:,:,aindex,kindex])

        plotting.plot_real_grid(rdata, title='Real part of wavelet coefficients', mask=mask, save=save, clabel=self.units, \
            contours=contours, **kwargs)
        plotting.plot_real_grid(idata, title='Imaginary part of wavelet coefficients', mask=mask, save=save, clabel=self.units, \
            contours=contours, **kwargs)


    def plot_scalogram(self, kindex=None, log=True, mask=None, title='Wavelet scalogram', \
        save=None, clabel=None, contours=None, **kwargs):
        """
        This method plots the wavelet scalogram of a :class:`~plateflex.classes.Grid` 
        object at a wavenumber index (int). Raises ``Exception`` for the cases where:

        - no wavenumber index is specified (kindex)
        - wavenumber index is lower than 0 or larger than self.ns - 1

        .. note::

            If no ``wl_sg`` attribute is found, the method automatically calculates
            the wavelet scalogram (and maybe also the wavelet transform) first.

        """

        if kindex is None:
            raise(Exception('Specify index of wavenumber for plotting'))

        if kindex>self.ns-1 or kindex<0:
            raise(Exception('Invalid index: should be between 0 and '+str(self.ns)))

        try:
            data = self.wl_sg[:,:,kindex]
        except:
            print('Calculating the scalogram first')
            self.wlet_scalogram()
            data = self.wl_sg[:,:,kindex]

        if log:
            plotting.plot_real_grid(data, log=log, mask=mask, title=title, save=save, \
                clabel=self.logsg_units, contours=contours, **kwargs)
        else:
            plotting.plot_real_grid(data, log=log, mask=mask, title=title, save=save, \
                clabel=self.sg_units, contours=contours, **kwargs)


class GravGrid(Grid):
    """
    Basic grid class of :mod:`~plateflex` for gravity data that inherits 
    from :class:`~plateflex.classes.Grid`

    .. rubric:: Additional Attributes

    ``units``: str
        Units of Gravity anomaly (':math:`mGal`')
    ``sg_units``: str
        Units of wavelet PSD (scalogram) (':math:`mGal^2/|k|`')
    ``logsg_units``: str
        Units of log of wavelet PSD (log(scalogram)) (':math:`log(mGal^2/|k|)`')
    ``title``: str
        Descriptor for Gravity data

    .. rubric:: Example

    >>> import numpy as np
    >>> from plateflex import Grid, GravGrid
    >>> nn = 200; dd = 10.
    >>> gravgrid = GravGrid(np.random.randn(nn, nn), dd, dd)
    >>> isinstance(gravgrid, Grid)
    True
    """

    def __init__(self, grid, dx, dy):

        Grid.__init__(self, grid, dx, dy)
        self.units = 'mGal'
        self.sg_units = r'mGal$^2$/|k|'
        self.logsg_units = r'log(mGal$^2$/|k|)'
        self.title = 'Gravity anomaly'

class BougGrid(GravGrid):
    """
    Basic grid class of :mod:`~plateflex` for Bouguer gravity data that inherits 
    from :class:`~plateflex.classes.GravGrid`

    .. rubric:: Additional Attributes

    ``title``: str
        Descriptor for Bouguer gravity data

    .. rubric: Example

    >>> import numpy as np
    >>> from plateflex import Grid, BougGrid, GravGrid
    >>> nn = 200; dd = 10.
    >>> bouggrid = BougGrid(np.random.randn(nn, nn), dd, dd)
    >>> isinstance(bouggrid, GravGrid)
    True
    >>> isinstance(bouggrid, Grid)
    True

    """

    def __init__(self, grid, dx, dy):

        GravGrid.__init__(self, grid, dx, dy)
        self.title = 'Bouguer anomaly'

class FairGrid(GravGrid):
    """
    Basic grid class of :mod:`~plateflex` for Free-air gravity data that inherits 
    from :class:`~plateflex.classes.GravGrid`

    .. rubric:: Additional Attributes

    ``title``: str
        Descriptor for Free-air gravity data

    .. rubric: Example

    >>> import numpy as np
    >>> from plateflex import Grid, FairGrid, GravGrid
    >>> nn = 200; dd = 10.
    >>> fairgrid = FairGrid(np.random.randn(nn, nn), dd, dd)
    >>> isinstance(fairgrid, GravGrid)
    True
    >>> isinstance(fairgrid, Grid)
    True
    """

    def __init__(self, grid, dx, dy):

        GravGrid.__init__(self, grid, dx, dy)
        self.title = 'Free-air anomaly'

class TopoGrid(Grid):
    """
    Basic grid class of :mod:`~plateflex` for Topography data that inherits 
    from :class:`~plateflex.classes.Grid`

    .. rubric:: Additional Attributes

    ``units`` : str
        Units of Topography (':math:`m`')
    ``sg_units`` : str
        Units of wavelet PSD (scalogram) (':math:`m^2/|k|`')
    ``logsg_units`` : str
        Units of log of wavelet PSD (log(scalogram)) (':math:`log(m^2/|k|)`')
    ``title``: str
        Descriptor for Topography data

    .. rubric: Example
    
    >>> import numpy as np
    >>> from plateflex import Grid, TopoGrid
    >>> nn = 200; dd = 10.
    >>> topogrid = TopoGrid(np.random.randn(nn, nn), dd, dd)
    >>> isinstance(topogrid, Grid)
    True

    .. note::

        Automatically converts grid values to 'meters' if standard deviation is lower 
        than 20
    """

    def __init__(self, grid, dx, dy):

        Grid.__init__(self, grid, dx, dy)
        self.units = 'm'
        self.sg_units = r'm$^2$/|k|'
        self.logsg_units = r'log(m$^2$/|k|)'
        self.title = 'Topography'

        if np.std(self.data) < 20.:
            self.data *= 1.e3
                   
        water_depth = grid.copy()
        water_depth[grid>0.] = 0.

        self.water_depth = -1.*water_depth

    def filter_water_depth(self, sigma=10, returned=False):
        from skimage.filters import gaussian

        water_depth = gaussian(self.data, sigma=sigma)
        water_depth[self.data>0.] = 0.

        self.water_depth = -1.*water_depth
        
        if returned:
            return water_depth

    def plot_water_depth(self, mask=None, title=None, save=None, clabel=None, contours=None, **kwargs):

        if title is not None:
            plotting.plot_real_grid(self.water_depth, title=title, mask=mask, save=save, \
                clabel=self.units, contours=contours, **kwargs)
        else:
            plotting.plot_real_grid(self.water_depth, title='Water depth', mask=mask, save=save, \
                clabel=self.units, contours=contours, **kwargs)


class RhocGrid(Grid):
    """
    Basic grid class of :mod:`~plateflex` for crustal density data that inherits 
    from :class:`~plateflex.classes.Grid`

    .. rubric:: Additional Attributes

    ``units`` : str
        Units of Density (':math:`kg/m^3`')
    ``title``: str
        Descriptor for Density data

    .. note::

        This class should only be used to specify the density of the
        crust at each cell location. Although the :class:`~plateflex.classes.Grid`
        methods are still available, they are not useful in this context.

    """

    def __init__(self, grid, dx, dy):

        Grid.__init__(self, grid, dx, dy)
        self.units = r'kg/m$^3$'
        self.sg_units = None
        self.logsg_units = None
        self.title = 'Crustal density'

        if np.std(self.data) < 10.:
            self.data *= 1.e3

class ZcGrid(Grid):
    """
    Basic grid class of :mod:`~plateflex` for crustal thickness data that inherits 
    from :class:`~plateflex.classes.Grid`

    .. rubric:: Additional Attributes

    ``units`` : str
        Units of Crustal thickness (':math:`m`')
    ``title``: str
        Descriptor for Crustal thickness data

    .. note::

        This class should only be used to specify the thickness of the
        crust at each cell location. Although the :class:`~plateflex.classes.Grid`
        methods are still available, they are not useful in this context.
         
    """

    def __init__(self, grid, dx, dy):

        Grid.__init__(self, grid, dx, dy)
        self.units = 'm'
        self.sg_units = None
        self.logsg_units = None
        self.title = 'Crustal thickness'

        if np.std(self.data) < 20.:
            self.data *= 1.e3

class Project(object):
    """
    Container for :class:`~plateflex.classes.Grid` objects, with
    methods to calculate the wavelet admittance and coherence and 
    estimate flexural model parameters as well as plot various results. 

    :type grids: list of :class:`~plateflex.classes.Grid`, optional
    :param grids: Initial list of PlateFlex :class:`~plateflex.classes.Grid`
        objects.

    .. rubric:: Default Attributes

    ``grids``: List
        List of :class:`~plateflex.classes.Grid` objects
    ``inverse`` : str
        Type of inversion to perform. By default the type is `'L2'` for
        non-linear least-squares. Options are: `'L2'` or `'bayes'`
    ``mask`` : Array
        2D array of boolean values determined independently
    ``initialized`` : Bool
        Whether or not the project has been initialized and is ready for 
        the wavelet analysis and estimation steps. By default this parameter
        is ``False``, unless the method :func:`~plateflex.classes.Project.init` 
        has been executed.

    .. note::

        Can hold a list of any length with any type of 
        :class:`~plateflex.classes.Grid` objects - however the wavelet 
        calculations will only proceed if the project holds one each 
        of :class:`~plateflex.classes.TopoGrid` and 
        :class:`~plateflex.classes.GravGrid` (``BougGrid`` or ``FairGrid``) objects. 

    .. rubric:: Examples

    Create zero-valued square grid

    >>> import numpy as np
    >>> from plateflex import Grid, TopoGrid, BougGrid, GravGrid, Project
    >>> nn = 200; dd = 10.
    >>> topogrid = TopoGrid(np.random.randn(nn, nn), dd, dd)
    >>> bouggrid = BougGrid(np.random.randn(nn, nn), dd, dd)
    >>> isinstance(bouggrid, GravGrid)
    True
    >>> isinstance(bouggrid, Grid)
    True

    Assign project with list of grids

    >>> Project(grids=[topogrid, bouggrid])
    <plateflex.classes.Project object at 0x10c62b320>

    Add Grid to project

    >>> project = Project()
    >>> project += topogrid
    >>> project += bouggrid
    >>> project.grids
    [<plateflex.classes.TopoGrid object at 0x1176b4240>, <plateflex.classes.BougGrid object at 0x1176b42e8>]

    Append Grid to project

    >>> project = Project()
    >>> project.append(topogrid)
    <plateflex.classes.Project object at 0x1176b9400>
    >>> project.grids[0]
    <plateflex.classes.TopoGrid object at 0x1176b4240>

    Initialize project

    >>> project.init()

    Exception if project does not contain exactly one TopoGrid and one GravGrid

    >>> project = Project(grids=[topogrid, topogrid])
    >>> project.init()
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
      File "/Users/pascalaudet/Softwares/Python/plateflex-dev/PlateFlex/plateflex/classes.py", line 492, in wlet_admit_coh
        grids = self.grids + other.grids
    Exception: There needs to be one GravGrid object in Project

    Calculate wavelet admittance and coherence

    >>> project = Project(grids=[topogrid, bouggrid])
    >>> project.wlet_admit_coh()
     #loops = 17:  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17
     #loops = 17:  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17
     Calculation jackknife error on admittance and coherence

    """

    def __init__(self, grids=None):

        self.inverse = 'L2'
        self.grids = []
        self.mask = None
        self.initialized = False

        if isinstance(grids, Grid):
            grids = [grids]
        if grids:
            self.grids.extend(grids)

    def __add__(self, other):
        """
        Add two `:class:`~plateflex.classes.Grid` objects or a 
        :class:`~plateflex.classes.Project` object with a single grid.

        """
        if isinstance(other, Grid):
            other = Project([other])
        if not isinstance(other, Project):
            raise TypeError
        grids = self.grids + other.grids
        return self.__class__(grids=grids)

    def __iter__(self):
        """
        Return a robust iterator for :class:`~plateflex.classes.Grid` 
        objects

        """
        return list(self.grids).__iter__()

    def append(self, grid):
        """
        Append a single :class:`~plateflex.classes.Grid` object to the 
        current `:class:`~plateflex.classes.Project` object.

        :type grid: :class:`~plateflex.classes.Grid`
        :param grid: object to append to project

        .. rubric:: Example
            
        >>> import numpy as np
        >>> from plateflex import Grid, Project
        >>> nn = 200; dd = 10.
        >>> grid = Grid(np.random.randn(nn, nn), dd, dd)
        >>> project = Project()
        >>> project.append(grid)
        """
        
        if isinstance(grid, Grid):
            self.grids.append(grid)
        else:
            msg = 'Append only supports a single Grid object as an argument.'
            raise TypeError(msg)

        return self

    def extend(self, grid_list):
        """
        Extend the current Project object with a list of Grid objects.

        :param trace_list: list of :class:`~plateflex.classes.Grid` objects or
            :class:`~plateflex.classes.Project`.

        .. rubric:: Example

        >>> import numpy as np
        >>> from plateflex import Grid, Project
        >>> nn = 200; dd = 10.
        >>> grid1 = Grid(np.random.randn(nn, nn), dd, dd)
        >>> grid2 = Grid(np.random.randn(nn, nn), dd, dd)
        >>> project = Project()
        >>> project.extend(grids=[grid1, grid2])

        """
        if isinstance(grid_list, list):
            for _i in grid_list:
                # Make sure each item in the list is a Grid object.
                if not isinstance(_i, Grid):
                    msg = 'Extend only accepts a list of Grid objects.'
                    raise TypeError(msg)
            self.grids.extend(grid_list)
        elif isinstance(grid_list, Project):
            self.grids.extend(grid_list.grids)
        else:
            msg = 'Extend only supports a list of Grid objects as argument.'
            raise TypeError(msg)
        return self


    def init(self):
        """
        Method to initialize a project. This step is required before calculating
        the wavelet admittance and coherence. The method checks that the project contains
        one each of :class:`~plateflex.classes.TopoGrid` and :class:`~plateflex.classes.GravGrid`.
        It also ensures that all grids have the same shape and sampling intervals. 
        If Grids of type :class:`~plateflex.classes.RhocGrid` and/or :class:`~plateflex.classes.ZcGrid`
        are present, the project attributes will be updated with data from those grids to be
        used in the estimation part. The method further sets the water depth at each cell location
        (grid point) from the :class:`~plateflex.classes.TopoGrid` object. 

        .. rubric:: Additional Attributes

        ``water_depth`` : :class:`~numpy.ndarray`
            Grid of water depth from topography data (shape (`nx,ny`))
        ``nx`` : int 
            Number of grid cells in the x-direction
        ``ny`` : int 
            Number of grid cells in the y-direction
        ``ns`` :  int 
            Number of wavenumber samples
        ``k`` : np.ndarray 
            1D array of wavenumbers
        ``initialized`` : bool
            Set to ``True`` when method is called successfully
            
        .. rubric:: Optional Attributes

        ``rhoc`` : :class:`~numpy.ndarray`
            Grid of crustal density data (shape (`nx,ny`))
        ``zc`` : :class:`~numpy.ndarray`
            Grid of crustal thickness data (shape (`nx,ny`))

        .. rubric:: Example

        >>> project = Project[grids=[topogrid, fairgrid, zcgrid]]
        >>> project.init()

        """
        

        # Methods will fail if there is no ``TopoGrid`` object in list
        if not any(isinstance(g, TopoGrid) for g in self.grids):
            raise(Exception('There needs to be one TopoGrid object in Project'))

        # Methods will fail if there is no ``GravGrid`` object in list
        if not any(isinstance(g, GravGrid) for g in self.grids):
            raise(Exception('There needs to be one GravGrid object in Project'))

        # Check that all grids have the same shape 
        shape = [grid.data.shape for grid in self.grids]
        if not (len(set(shape))==1):
            raise(Exception('Grids do not have the same shape - aborting:'+str(shape)))

        # Check that all grids have the same sampling intervals
        dd = [(grid.dx, grid.dy) for grid in self.grids]
        if not (len(set(dd))==1):
            raise(Exception('Grids do not have the same sampling intervals - aborting:'+str(dd)))

        # Check whether gravity grid is Free air or Bouguer and set global variable accordingly
        if any(isinstance(g, BougGrid) for g in self.grids):
            cf_f.boug = 1
        elif any(isinstance(g, FairGrid) for g in self.grids):
            cf_f.boug = 0

        # Initialize model attributes to None (i.e., default values will be used)
        self.rhoc = None
        self.zc = None

        # Identify the ``Grid`` types and set new attributes if available
        for grid in self.grids:
            if isinstance(grid, RhocGrid):
                self.rhoc = grid.data
            elif isinstance(grid, ZcGrid):
                self.zc = grid.data
            elif isinstance(grid, TopoGrid):
                self.water_depth = grid.water_depth

        # Now set project attributes from first grid object
        self.k = self.grids[0].k
        self.ns = self.grids[0].ns
        self.nx = self.grids[0].nx
        self.ny = self.grids[0].ny
        self.initialized = True


    def wlet_admit_coh(self):
        """
        This method uses the module :mod:`~plateflex.cpwt.cpwt` to calculate 
        the wavelet admittance and coherence. The object needs to contain 
        exactly two :class:`~plateflex.classes.Grid` objects, one of each of 
        :class:`~plateflex.classes.TopoGrid` and :class:`~plateflex.classes.GravGrid` 
        objects. 

        .. note::

            If no ``wl_trans`` attribute is found for individual 
            :class:`~plateflex.classes.Grid` objects, the method automatically calculates
            them first.

        Stores the wavelet admittance, coherence and their errors as attributes

        .. rubric:: Additional Attributes

        ``wl_admit`` : :class:`~numpy.ndarray`
            Wavelet admittance (shape (`nx,ny,ns`))
        ``wl_eadmit`` : :class:`~numpy.ndarray`
            Error of wavelet admittance (shape (`nx,ny,ns`))
        ``wl_coh`` : :class:`~numpy.ndarray`
            Wavelet coherence (shape (`nx,ny,ns`))
        ``wl_ecoh`` : :class:`~numpy.ndarray`
            Error of wavelet coherence (shape (`nx,ny,ns`))

        """

        if not self.initialized:
            raise(Exception("Project not yet initialized - Abort"))

        # Identify the ``Grid`` types for proper calculation of admittance and coherence 
        for grid in self.grids:
            if isinstance(grid, TopoGrid):
                try:
                    wl_trans_topo = grid.wl_trans
                except:
                    grid.wlet_transform()
                    wl_trans_topo = grid.wl_trans
            elif isinstance(grid, GravGrid):
                try:
                    wl_trans_grav = grid.wl_trans
                except:
                    grid.wlet_transform()
                    wl_trans_grav = grid.wl_trans

        # Calculate wavelet admittance and coherence by calling function wlet_admit_coh
        # from module ``plateflex.cpwt.cpwt``
        wl_admit, ewl_admit, wl_coh, ewl_coh = \
            cpwt.wlet_admit_coh(wl_trans_topo, wl_trans_grav)

        # Store the admittance, coherence and error grids as attributes
        self.wl_admit = wl_admit
        self.ewl_admit = ewl_admit
        self.wl_coh = wl_coh
        self.ewl_coh = ewl_coh

        return 

    def plot_admit_coh(self, kindex=None, mask=None, save=None, contours=None, **kwargs):
        """
        Method to plot grids of wavelet admittance and coherence at a given  wavenumber index. 

        :type kindex: int
        :param kindex: Index of wavenumber array
        :type mask: :class:`~numpy.ndarray`, optional 
        :param mask: Array of booleans for masking data points
        :type title: str, optional
        :param title: Title of plot
        :type save: str, optional
        :param save: Name of file for to save figure
        :type clabel: str, optional    
        :param clabel: Label for colorbar
        :type contours: List
        :param contours: List of contour lines obtained separately
        :type kwargs: Keyword arguments
        :param kwargs: Keyword arguments allowing more control on plots
        """

        if kindex is None:
            raise(Exception('Specify index of wavenumber for plotting'))

        try:
            adm = self.wl_admit[:,:,kindex]
            coh = self.wl_coh[:,:,kindex]
        except:
            print('Calculating the admittance and coherence first')
            self.wlet_admit_coh()
            adm = self.wl_admit[:,:,kindex]
            coh = self.wl_coh[:,:,kindex]

        if kindex>self.ns or kindex<0:
            raise(Exception('Invalid index: should be between 0 and '+str(self.ns)))

        plotting.plot_real_grid(adm, mask=mask, title='Admittance', save=save, clabel='mGal/m', \
            contours=contours, **kwargs)
        plotting.plot_real_grid(coh, mask=mask, title='Coherence', save=save, clabel=None, \
            contours=contours, **kwargs)


    def estimate_cell(self, cell=(0,0), alph=False, atype='joint', returned=False):
        """
        Method to estimate the parameters of the flexural model at a single cell location
        of the input grids. The type of estimation performed is set by the project attribute 
        ``inverse``. See :class:`~plateflex.classes.Project`` and :mod:`~plateflex.estimate` 
        for details. The model parameters ``rhoc`` and ``zc`` are extracted from the 
        corresponding :class:`~plateflex.classes.RhocGrid` and :class:`~plateflex.classes.ZcGrid`
        objects if they were initialized in the project.

        :type cell: tuple
        :param cell: Indices of cell location within grid
        :type alph: bool, optional
        :param alph: Whether or not to estimate parameter ``alpha``
        :type atype: str, optional
        :param atype: Whether to use the admittance ('admit'), coherence ('coh') or both ('joint')
        :type returned: bool, optional
        :param returned: Whether or not to return the estimates

        .. rubric:: Additional Attributes

        ``trace`` : :class:`~pymc3.backends.base.MultiTrace`
            Posterior samples from the MCMC chains
        ``summary`` : :class:`~pandas.core.frame.DataFrame`
            Summary statistics from Posterior distributions
        ``map_estimate`` : dict, optional
            Container for Maximum a Posteriori (MAP) estimates
        ``cell`` : tuple 
            Indices of cell location within grid
        
        Results are stored as attributes of :class:`~plateflex.classes.Project` 
        object.
        """

        if not isinstance(alph, bool):
            raise(Exception("'alph' should be a boolean: defaults to False"))

        if atype not in ['admit', 'coh', 'joint']:
            raise(Exception("'atype' should be one among: 'admit', 'coh', or 'joint'"))

        # Extract admittance and coherence at cell indices
        adm = self.wl_admit[cell[0], cell[1], :]
        eadm = self.ewl_admit[cell[0], cell[1], :]
        coh = self.wl_coh[cell[0], cell[1], :]
        ecoh = self.ewl_coh[cell[0], cell[1], :]

        # Now set model parameters if they are availble as attributes
        cf_f.wd = self.water_depth[cell[0], cell[1]]
        if self.rhoc is not None:
            cf_f.rhoc = self.rhoc[cell[0], cell[1]]
        if self.zc is not None:
            cf_f.zc = self.zc[cell[0], cell[1]]

        if self.inverse=='L2':
            summary = estimate.L2_estimate_cell( \
                self.k, adm, eadm, coh, ecoh, alph, atype)

            # Return estimates if requested
            if returned:
                return summary

            # Otherwise store as object attributes
            else:
                self.alph = alph
                self.atype = atype
                self.cell = cell
                self.summary = summary

        elif self.inverse=='bayes':
            trace, summary, map_estimate = estimate.bayes_estimate_cell( \
                self.k, adm, eadm, coh, ecoh, alph, atype)

            # Return estimates if requested
            if returned:
                return summary, map_estimate

            # Otherwise store as object attributes
            else:
                self.alph = alph
                self.atype = atype
                self.cell = cell
                self.trace = trace
                self.map_estimate = map_estimate
                self.summary = summary

    def estimate_grid(self, nn=10, alph=False, atype='joint', parallel=False):
        """
        Method to estimate the parameters of the flexural model at all grid point locations.
        It is also possible to decimate the number of grid cells at which to estimate parameters. 

        :type nn: int
        :param nn: Decimator. If grid shape is ``(nx, ny)``, resulting grids will have 
            shape of ``(int(nx/nn), int(ny/nn))``. 
        :type alph: bool, optional 
        :param alph: Whether or not to estimate parameter ``alpha``
        :type atype: str, optional
        :param atype: Whether to use the admittance (`'admit'`), coherence (`'coh'`) or 
            both (`'joint'`)

        .. rubric:: Additional Attributes

        ``mean_Te_grid`` : :class:`~numpy.ndarray` 
            Grid with mean Te estimates (shape ``(nx, ny``))
        ``MAP_Te_grid`` : :class:`~numpy.ndarray` 
            Grid with MAP Te estimates (shape ``(nx, ny``))
        ``std_Te_grid`` : :class:`~numpy.ndarray` 
            Grid with std Te estimates (shape ``(nx, ny``))
        ``mean_F_grid`` : :class:`~numpy.ndarray` 
            Grid with mean F estimates (shape ``(nx, ny``))
        ``MAP_F_grid`` : :class:`~numpy.ndarray` 
            Grid with MAP F estimates (shape ``(nx, ny``))
        ``std_F_grid`` : :class:`~numpy.ndarray` 
            Grid with std F estimates (shape ``(nx, ny``))

        .. rubric:: Optional Attributes

        ``mean_a_grid`` : :class:`~numpy.ndarray` 
            Grid with mean alpha estimates (shape ``(nx, ny``)). Only present if ``alph=True``.
        ``MAP_a_grid`` : :class:`~numpy.ndarray` 
            Grid with MAP alpha estimates (shape ``(nx, ny``)). Only present if ``alph=True``.
        ``std_a_grid`` : :class:`~numpy.ndarray` 
            Grid with std alpha estimates (shape ``(nx, ny``)). Only present if ``alph=True``.
        ``chi2_grid`` : :class:`~numpy.ndarray` 
            Grid with reduced chi-squared estimates (shape ``(nx, ny``)). 
            Only present if ``project.inverse='L2'``
        ``new_mask`` : Array
            New grid of masked (boolean) values, corresponding to decimated mask.

        """

        self.nn = nn

        if not isinstance(alph, bool):
            raise(Exception("'alph' should be a boolean: defaults to False"))
        self.alph = alph

        if atype not in ['admit', 'coh', 'joint']:
            raise(Exception("'atype' should be one among: 'admit', 'coh', or 'joint'"))
        self.atype = atype

        # Initialize result grids to zeroes
        if self.mask is not None:
            new_mask_grid = np.zeros((int(self.nx/nn),int(self.ny/nn)), dtype=bool)

        mean_Te_grid = np.zeros((int(self.nx/nn),int(self.ny/nn)))
        std_Te_grid = np.zeros((int(self.nx/nn),int(self.ny/nn)))
        mean_F_grid = np.zeros((int(self.nx/nn),int(self.ny/nn)))
        std_F_grid = np.zeros((int(self.nx/nn),int(self.ny/nn)))
        if self.alph:
            mean_a_grid = np.zeros((int(self.nx/nn),int(self.ny/nn)))
            std_a_grid = np.zeros((int(self.nx/nn),int(self.ny/nn)))

        if self.inverse=='bayes':
            MAP_Te_grid = np.zeros((int(self.nx/nn),int(self.ny/nn)))
            MAP_F_grid = np.zeros((int(self.nx/nn),int(self.ny/nn)))
            if self.alph:
                MAP_a_grid = np.zeros((int(self.nx/nn),int(self.ny/nn)))

        elif self.inverse=='L2':
            chi2_grid = np.zeros((int(self.nx/nn),int(self.ny/nn)))

        if parallel:

            raise(Exception('Parallel implementation does not work - check again later'))
            # from joblib import Parallel, delayed
            # cf.cores=1

            # # Run nested for loop in parallel to cover the whole grid
            # results = Parallel(n_jobs=4)(delayed(self.estimate_cell) \
            #     (cell=(i,j), alph=alph, atype=atype, returned=True) \
            #     for i in range(0, self.nx-nn, nn) for j in range(0, self.ny-nn, nn))

        else:
            for i in _progressbar(range(0, self.nx-nn, nn), 'Computing: ', 10):
                for j in range(0, self.ny-nn, nn):
                    
                    # # For reference - index values
                    #print(i,j)
                    #print(self.mask[i,j])

                    # tuple of cell indices
                    cell = (i,j)

                    # Skip masked cells
                    if self.mask is not None:
                        new_mask_grid[int(i/nn),int(j/nn)] = self.mask[i,j]
                        if self.mask[i,j]:
                            continue

                    if self.inverse=='bayes':

                        # Carry out calculations by calling the ``estimate_cell`` method
                        summary, map_estimate = self.estimate_cell(cell=cell, \
                            alph=alph, atype=atype, returned=True)

                        # Extract estimates from summary and map_estimate
                        res = estimate.get_bayes_estimates(summary, map_estimate)

                        # Distribute the parameters back to space
                        mean_Te = res[0]
                        std_Te = res[1]
                        MAP_Te = res[4]
                        mean_F = res[5]
                        std_F = res[6]
                        MAP_F = res[9]
                        if self.alph:
                            mean_a = res[10]
                            std_a = res[11]
                            MAP_a = res[14]

                        # Store values in smaller arrays
                        mean_Te_grid[int(i/nn),int(j/nn)] = mean_Te
                        MAP_Te_grid[int(i/nn),int(j/nn)] = MAP_Te
                        std_Te_grid[int(i/nn),int(j/nn)] = std_Te
                        mean_F_grid[int(i/nn),int(j/nn)] = mean_F
                        MAP_F_grid[int(i/nn),int(j/nn)] = MAP_F
                        std_F_grid[int(i/nn),int(j/nn)] = std_F
                        if self.alph:
                            mean_a_grid[int(i/nn),int(j/nn)] = mean_a
                            MAP_a_grid[int(i/nn),int(j/nn)] = MAP_a
                            std_a_grid[int(i/nn),int(j/nn)] = std_a

                    elif self.inverse=='L2':
                        #print('l2')

                        # Carry out calculations by calling the ``estimate_cell`` method
                        summary = self.estimate_cell(cell=cell, \
                            alph=alph, atype=atype, returned=True)
                        print(summary)

                        # Extract estimates from summary and map_estimate
                        res = estimate.get_L2_estimates(summary)

                        # Distribute the parameters back to space
                        mean_Te = res[0]
                        std_Te = res[1]
                        mean_F = res[2]
                        std_F = res[3]
                        if self.alph:
                            mean_a = res[4]
                            std_a = res[5]
                            chi2 = res[6]
                        else:
                            chi2 = res[4]

                        # Store values in smaller arrays
                        mean_Te_grid[int(i/nn),int(j/nn)] = mean_Te
                        std_Te_grid[int(i/nn),int(j/nn)] = std_Te
                        mean_F_grid[int(i/nn),int(j/nn)] = mean_F
                        std_F_grid[int(i/nn),int(j/nn)] = std_F
                        chi2_grid[int(i/nn),int(j/nn)] = chi2
                        if self.alph:
                            mean_a_grid[int(i/nn),int(j/nn)] = mean_a
                            std_a_grid[int(i/nn),int(j/nn)] = std_a

        if self.mask is not None:
            self.new_mask_grid = new_mask_grid

        if self.inverse=='bayes':

            # Store grids as attributes
            self.mean_Te_grid = mean_Te_grid
            self.MAP_Te_grid = MAP_Te_grid
            self.std_Te_grid = std_Te_grid
            self.mean_F_grid = mean_F_grid
            self.MAP_F_grid = MAP_F_grid
            self.std_F_grid = std_F_grid
            if self.alph:
                self.mean_a_grid = mean_a_grid
                self.MAP_a_grid = MAP_a_grid
                self.std_a_grid = std_a_grid

        elif self.inverse=='L2':
            # Store grids as attributes
            self.chi2_grid = chi2_grid
            self.mean_Te_grid = mean_Te_grid
            self.std_Te_grid = std_Te_grid
            self.mean_F_grid = mean_F_grid
            self.std_F_grid = std_F_grid
            if self.alph:
                self.mean_a_grid = mean_a_grid
                self.std_a_grid = std_a_grid


    def plot_bayes_stats(self, title=None, save=None):
        """
        Method to plot the marginal and joint distributions of samples drawn from the 
        posterior distribution as well as the extracted statistics. Calls the function 
        :func:`~plateflex.plotting.plot_stats` with attributes as arguments.

        :type title: str, optional 
        :param title: Title of plot
        :type save: str, optional
        :param save: Name of file for to save figure

        """

        try:
            plotting.plot_bayes_stats(self.trace, self.summary, \
                self.map_estimate, title=title, save=save)
        except:
            raise(Exception("No 'cell' estimate available"))

    def plot_functions(self, est='mean', title=None, save=None):
        """
        Method to plot observed and fitted admittance and coherence functions using 
        one of ``mean`` or ``MAP`` estimates. The ``MAP`` is only available if the project
        attribute has been set to ``project='bayes'``. Calls the function 
        :func:`~plateflex.plotting.plot_functions` with attributes as arguments.

        :type est: str, optional
        :param est: Type of inference estimate to use for predicting admittance and coherence
        :type title: str, optional 
        :param title: Title of plot
        :type save: str, optional
        :param save: Name of file for to save figure
        """

        if est not in ['mean', 'MAP']:
            raise(Exception("Choose one among: 'mean', or 'MAP'"))
            
        try:
            cell = self.cell
            k = self.k
            adm = self.wl_admit[cell[0], cell[1], :]
            eadm = self.ewl_admit[cell[0], cell[1], :]
            coh = self.wl_coh[cell[0], cell[1], :]
            ecoh = self.ewl_coh[cell[0], cell[1], :]
            cf_f.wd = self.water_depth[cell[0], cell[1]]

        except:
            raise(Exception("No cell specified. Aborting"))

        try:
            ma = np.pi/2.

            if self.inverse=='bayes':
                # Extract statistics from summary object
                if est=='mean':
                    mte = self.summary.loc['Te', est]
                    mF = self.summary.loc['F', est]
                    if sum(self.summary.index.isin(['alpha']))==1:
                        ma = self.summary.loc['alpha', est]

                # Extract MAP from map_estimate object
                elif est=='MAP':
                    mte = np.float(self.map_estimate['Te'])
                    mF = np.float(self.map_estimate['F'])
                    if 'alpha' in self.map_estimate:
                        ma = np.float(self.map_estimate['alpha'])
                else:
                    raise(Exception("estimate does not exist. Choose among: 'mean' or 'MAP'"))

            elif self.inverse=='L2':

                # Extract statistics from summary object
                mte = self.summary.loc['Te', 'mean']
                mF = self.summary.loc['F', 'mean']
                if sum(self.summary.index.isin(['alpha']))==1:
                    ma = self.summary.loc['alpha', 'mean']

            # Calculate predicted admittance and coherence from estimates
            padm, pcoh = estimate.real_xspec_functions(k, mte, mF, ma)

            # Call function from ``plotting`` module
            plotting.plot_functions(k, adm, eadm, coh, ecoh, \
                padm=padm, pcoh=pcoh, title=title, save=save)

        except:

            # Call function from ``plotting`` module
            plotting.plot_functions(k, adm, eadm, coh, ecoh, \
                title=title, save=save)


    def plot_results(self, mean_Te=False, MAP_Te=False, std_Te=False, \
        mean_F=False, MAP_F=False, std_F=False, mean_a=False, MAP_a=False, \
        std_a=False, chi2=False, mask=False, contours=None, save=None, \
        filter=True, sigma=1, **kwargs):
        """
        Method to plot grids of estimated parameters with fixed labels and titles. 
        To have more control over the plot rendering, use the function 
        :func:`~plateflex.plotting.plot_real_grid` with the relevant quantities and 
        plotting options.

        :type mean/MAP/std_Te/F/a: bool
        :param mean/MAP/std_Te/F/a: Type of plot to produce. 
            All variables default to False (no plot generated)
        :type mask: bool
        :param mask: Whether or not to plot the mask
        :type contours: List
        :param contours: List of contours with coordinate positions
        :type filter: bool
        :param filter: Whether or not to filter the resulting grid using a Gaussian filter
        :type sigma: int
        :param sigma: Standard deviation of filter (in terms of adjacent cells), if set to ``True``
        :type kwargs: Keyword arguments
        :param kwargs: Keyword arguments allowing more control on plots

        .. note::

           It is advisable to plot each grid separately using one call per grid, in order
           to have more control on plot paramters (e.g., colormap, min and max values of 
           colorbar, etc.).

        """

        from skimage.filters import gaussian

        if mask:
            try:
                new_mask = self.new_mask_grid
            except:
                new_mask = None
                print('No new mask found. Plotting without mask')
        else:
            new_mask=None

        if contours is not None:
            contours = np.array(contours)/self.nn

        if mean_Te:
            if filter:
                mean_Te_grid = gaussian(self.mean_Te_grid, sigma=sigma)
            else:
                mean_Te_grid = self.mean_Te_grid
            plotting.plot_real_grid(mean_Te_grid, mask=new_mask, \
                title='Elastic thickness', clabel='$T_e$ (km)', contours=contours, \
                save=save, **kwargs)
        if MAP_Te:
            if filter:
                MAP_Te_grid = gaussian(self.MAP_Te_grid, sigma=sigma)
            else:
                MAP_Te_grid = self.MAP_Te_grid
            plotting.plot_real_grid(MAP_Te_grid, mask=new_mask, \
                title='MAP estimate of $T_e$', clabel='$T_e$ (km)', contours=contours, \
                save=save, **kwargs)
        if std_Te:
            if filter:
                std_Te_grid = gaussian(self.std_Te_grid, sigma=sigma)
            else:
                std_Te_grid = self.std_Te_grid
            plotting.plot_real_grid(std_Te_grid, mask=new_mask, \
                title='Error on $T_e$', clabel='$T_e$ (km)', contours=contours, \
                save=save, **kwargs)
        if mean_F:
            if filter:
                mean_F_grid = gaussian(self.mean_F_grid, sigma=sigma)
            else:
                mean_F_grid = self.mean_F_grid
            plotting.plot_real_grid(mean_F_grid, mask=new_mask, \
                title='Initial load ratio', clabel='$F$', contours=contours, \
                save=save, **kwargs)
        if MAP_F:
            if filter:
                MAP_F_grid = gaussian(self.MAP_F_grid, sigma=sigma)
            else:
                MAP_F_grid = self.MAP_F_grid
            plotting.plot_real_grid(MAP_F_grid, mask=new_mask, \
                title='MAP estimate of $F$', clabel='$F$', contours=contours, \
                save=save, **kwargs)
        if std_F:
            if filter:
                std_F_grid = gaussian(self.std_F_grid, sigma=sigma)
            else:
                std_F_grid = self.std_F_grid       
            plotting.plot_real_grid(std_F_grid, mask=new_mask, \
                title='Error on $F$', clabel='$F$', contours=contours, \
                save=save, **kwargs)
        if mean_a:
            try:
                if filter:
                    mean_a_grid = gaussian(self.mean_a_grid, sigma)
                else:
                    mean_a_grid = self.mean_a_grid
                plotting.plot_real_grid(mean_a_grid, mask=new_mask, \
                    title=r'Phase difference between initial loads', clabel=r'$\alpha$', \
                    contours=contours, save=save, **kwargs)
            except:
                print("parameter 'alpha' was not estimated")
        if MAP_a:
            try:
                if filter:
                    MAP_a_grid = gaussian(self.MAP_a_grid, sigma=sigma)
                else:
                    MAP_a_grid = self.MAP_a_grid
                plotting.plot_real_grid(MAP_a_grid, mask=new_mask, \
                    title=r'MAP estimate of $\alpha$', clabel=r'$\alpha$', \
                    contours=contours, save=save, **kwargs)
            except:
                print("parameter 'alpha' was not estimated")
        if std_a:
            try:
                if filter:
                    std_a_grid = gaussian(self.std_a_grid, sigma=sigma)
                else:
                    std_a_grid = self.std_a_grid
                plotting.plot_real_grid(std_a_grid, mask=new_mask, \
                    title=r'Error on $\alpha$', clabel=r'$\alpha$', contours=contours, \
                    save=save, **kwargs)
            except:
                print("parameter 'alpha' was not estimated")
        if chi2:
            try:
                if filter:
                    chi2_grid = gaussian(self.chi2_grid, sigma=sigma)
                else:
                    chi2_grid = self.chi2_grid
                plotting.plot_real_grid(chi2_grid, mask=new_mask, \
                    title='Reduced chi-squared error', clabel=r'$\chi_{\nu}^2$', \
                    contours=contours, save=save, **kwargs)
            except:
                print("parameter 'chi2' was not estimated")



def _lam2k(nx, ny, dx, dy, p=0.85):
    """
    Calculate the optimal set of equally-spaced equivalent wavenumbers for given grid 
    parameters to be used in the wavelet analysis through the :mod:`~plateflex.cpwt.cpwt` 
    module. 

    :type nx: int
    :param nx: Size of grid in x direction
    :type ny: int
    :param ny: Size of grid in y direction
    :type dx: float 
    :param dx: Sample distance in x direction (km)
    :type dy: float 
    :param dy: Sample distance in y direction (km)
    :type p: float 
    :param p: Fractional amplitude at which wavelets overlap in wavenumber space

    :return:  
        (tuple): Tuple containing:
            * ns (int): Size of wavenumber array
            * k (:class:`~numpy.ndarray`): Wavenumbers (rad/m)

    .. note::

        This is different from the exact Fourier wavenumbers calculated for 
        the grid using :func:`~numpy.fft.fftfreq`, as the continuous wavelet
        transform can be defined at arbitrary wavenumbers.

    .. rubric:: Example

    >>> from plateflex import classes
    >>> # Define fake grid
    >>> nx = ny = 300
    >>> dx = dy = 20.
    >>> classes._lam2k(nx, ny, dx, dy)
    (18, array([2.96192196e-06, 3.67056506e-06, 4.54875180e-06, 5.63704569e-06,
           6.98571509e-06, 8.65705514e-06, 1.07282651e-05, 1.32950144e-05,
           1.64758613e-05, 2.04177293e-05, 2.53026936e-05, 3.13563910e-05,
           3.88584421e-05, 4.81553673e-05, 5.96765921e-05, 7.39542826e-05,
           9.16479263e-05, 1.13574794e-04]))
    """

    # Max is quarter of maximum possible wavelength from grid size
    maxlam = np.sqrt((nx*dx)**2. + (ny*dy)**2.)/4.*1.e3
    
    # Min is twice the minimum possible wavelength from grid size
    minlam = np.sqrt((2.*dx)**2. + (2.*dy)**2.)*1.e3

    # Assign first wavenumber
    lam = []; k = []; s = []
    lam.append(maxlam)
    k.append(2.*np.pi/lam[0])
    s.append(cf_w.k0/k[0])
    dk = np.sqrt(-2.*np.log(p))/s[0]
    ss = 1

    # Loop through k's until minlam is reached
    while lam[ss-1]>minlam:
            s = (cf_w.k0-np.sqrt(-2.*np.log(p)))/(k[ss-1]+dk)
            k.append(cf_w.k0/s)
            lam.append(2.*np.pi/k[ss])
            dk = np.sqrt(-2.0*np.log(p))/s
            ss = ss + 1
    ns = ss

    # Compute wavenumbers
    lam = np.array(lam)
    k = np.array(2.*np.pi/lam)

    return ns, k

import sys

def _progressbar(it, prefix="", size=60, file=sys.stdout):
    count = len(it)
    def show(j):
        x = int(size*j/count)
        file.write("%s[%s%s] %i/%i\r" % (prefix, "#"*x, "."*(size-x), j, count))
        file.flush()        
    show(0)
    for i, item in enumerate(it):
        yield item
        show(i+1)
    file.write("\n")
    file.flush()