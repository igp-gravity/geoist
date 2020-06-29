"""
Misc utilities
"""
from pathlib import Path
import sys
import hashlib
import h5py
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import ImageGrid
from mpl_toolkits.basemap import Basemap
import pyproj
from typing import Tuple, Any
import numpy as np
import datetime

import matplotlib.pyplot as plt
#from win32gui import GetWindowText, GetForegroundWindow
from PIL import Image
#import win32clipboard
import io

from datetime import timedelta
import xarray
from matplotlib.dates import DateFormatter, MinuteLocator, SecondLocator

def os_cache(project, platform=None):
    """
    Default cache location based on the operating system.

    Will insert the project name in the proper location of the path.

    Parameters
    ----------
    project : str
        The project name.
    platform : str or None
        The name of operating system as returned by ``sys.platform`` (``'darwin'`` for
        Mac, ``'win32'`` for Windows, and anything else will be treated as generic
        Linux/Unix. If None, will use the value of ``sys.platform``.

    Returns
    -------
    cache_path : :class:`pathlib.Path`
        The default location for the data cache. User directories (``'~'``) are not
        expanded.

    Examples
    --------

    >>> for os in ['darwin', 'win32', 'anything else']:
    ...     path = os_cache("myproject", platform=os)
    ...     print(path.parts)
    ('~', 'Library', 'Caches', 'myproject')
    ('~', 'AppData', 'Local', 'myproject', 'cache')
    ('~', '.cache', 'myproject')

    """
    if platform is None:
        platform = sys.platform
    if platform == "darwin":
        cache_path = Path("~", "Library", "Caches", project)
    elif platform == "win32":
        cache_path = Path("~", "AppData", "Local", project, "cache")
    else:  # *NIX
        cache_path = Path("~", ".cache", project)
    return cache_path


def file_hash(fname):
    """
    Calculate the SHA256 hash of a given file.

    Useful for checking if a file has changed or been corrupted.

    Parameters
    ----------
    fname : str
        The name of the file.

    Returns
    -------
    hash : str
        The hash of the file.

    Examples
    --------

    >>> fname = "test-file-for-hash.txt"
    >>> with open(fname, "w") as f:
    ...     __ = f.write("content of the file")
    >>> print(file_hash(fname))
    0fc74468e6a9a829f103d069aeb2bb4f8646bad58bf146bb0e3379b759ec4a00
    >>> import os
    >>> os.remove(fname)

    """
    # Calculate the hash in chunks to avoid overloading the memory
    chunksize = 65536
    hasher = hashlib.sha256()
    with open(fname, "rb") as fin:
        buff = fin.read(chunksize)
        while buff:
            hasher.update(buff)
            buff = fin.read(chunksize)
    return hasher.hexdigest()

def make_registry(directory, output, recursive=True):
    """
    Make a registry of files and hashes for the given directory.

    This is helpful if you have many files in your test dataset as it keeps you
    from needing to manually update the registry.

    Parameters
    ----------
    directory : str
        Directory of the test data to put in the registry. All file names in the
        registry will be relative to this directory.
    output : str
        Name of the output registry file.
    recursive : bool
        If True, will recursively look for files in subdirectories of *directory*.

    """
    directory = Path(directory)
    if recursive:
        pattern = "**/*"
    else:
        pattern = "*"

    files = sorted(
        [
            str(path.relative_to(directory))
            for path in directory.glob(pattern)
            if path.is_file()
        ]
    )

    hashes = [file_hash(str(directory / fname)) for fname in files]

    with open(output, "w") as outfile:
        for fname, fhash in zip(files, hashes):
            # Only use Unix separators for the registry so that we don't go insane
            # dealing with file paths.
            outfile.write("{} {}\n".format(fname.replace("\\", "/"), fhash))

def csr2hdf5(matrix,filename,mode='w',group_name='Mcsr'):
    with h5py.File(filename,mode) as f:
        g = f.create_group(group_name)
        g.create_dataset('data',data=matrix.data)
        g.create_dataset('indptr',data=matrix.indptr)
        g.create_dataset('indices',data=matrix.indices)
        g.attrs['shape'] = matrix.shape

def hdf52csr(filename,mode='r',group_name='Mcsr'):
    with h5py.File(filename,mode) as f:
        g = f[group_name]
        M1 = csr_matrix((g['data'][:],g['indices'][:],g['indptr'][:]),
                        g.attrs['shape'])
    return(M1)

def get_prism_pos(index,nx,ny,nz):
    '''calculate prism's position given its index.

    Args:
        index (int): prism's index
        nx,ny,nz (int): how many prisms along x,y,z axis.

    Returns:
        position of the prism
    '''
    k = index // (nx * ny)
    j = (index - k*(nx*ny))//nx
    i = (index - k*(nx*ny) - j*nx)
    return i,j,k

def get_prism_index(ix,iy,iz,nx,ny,nz):
    '''calculate prism's index given its position

    Args:
        ix,iy,iz (int): prism's position
        nx,ny,nz (int): how many prisms along x,y,z axis.

    Returns:
        int: index of the prism
    '''
    return iz*nx*ny + iy*nx + ix

def plot_kernel(ggz,nxyz=None,nobs=None,obs_extent=(-100,100,-100,100),image_grid=(3,5),fname=None):
    '''inspect the kernel matrix

    Args:
        ggz (ndarray): Kernel matrix. Each column correspond to a source point. Each row correspond to an observe station.
        nxyz (tuple): How many source points along x,y,z axis respectively.
        nobs (tuple): How many observe stations along x,y axis respectively.
        obs_extent (tuple): Define the observe area in a order of (min_x,max_x,min_y,max_y).
        image_grid (tuple): Define the dimension of image grid.
    '''
    if nxyz is None:
        nx = ny = nz = int(round(np.power(ggz.shape[1],1./3.)))
    else:
        nx,ny,nz = nxyz
    if nobs is None:
        obs_x = obs_y = int(round(np.power(ggz.shape[0],1./2.)))
    else:
        obs_x,obs_y = nobs
    n_rows,n_cols = image_grid
    fmt = ticker.ScalarFormatter()
    fmt.set_powerlimits((0,0))
    fig = plt.figure(figsize=(15,12))
    iz = np.linspace(0,nz-1,n_rows).astype(np.int32)
    ind = []
    for i in iz:
        ind.extend([i*nx*ny,i*nx*ny+nx-1,i*nx*ny+nx*ny//2+nx//2,(i+1)*nx*ny-nx,(i+1)*nx*ny-1])
    axes = ImageGrid(fig, 111,  # similar to subplot(111)
                     nrows_ncols=(n_rows,n_cols),
                     axes_pad=0.5,
                     add_all=True,
                     label_mode="L",
                     cbar_mode = 'each',
                     cbar_location = 'right',
                     cbar_pad='30%'
                     )
    i = 0
    for row in range(n_rows):
        for col in range(n_cols):
            ixyz = get_prism_pos(int(ind[i]),nx,ny,nz)
            im = axes[col+row*n_cols].imshow(ggz[:,int(ind[i])].reshape(-1,obs_x).transpose(),extent=obs_extent,origin='lower')
            axes[col+row*n_cols].set_title('layer {} of {}'.format(ixyz[2],nz))
            i += 1
            axes.cbar_axes[col+row*n_cols].colorbar(im,format=fmt)
    if fname is None:
        plt.show()
    else:
        plt.savefig(fname)

def plot_matrix(A,cbar_location='right',figsize=(18,18),cmap='coolwarm',fname=None):
    fig = plt.figure(figsize=figsize)
    axes = ImageGrid(fig, 111,  # similar to subplot(111)
                 nrows_ncols=(1,1),
                 axes_pad=2.0,
                 add_all=True,
                 label_mode="L",
                 cbar_mode = 'each',
                 cbar_location = cbar_location,
                 cbar_pad='2%'
                 )
    im = axes[0].imshow(A,cmap=cmap,interpolation='none')
    axes.cbar_axes[0].colorbar(im)
    if fname is None:
        plt.show()
    else:
        plt.savefig(fname)
    return fig,axes

def grid2Grid(x, y, data):
    from geoist.others.geodict import GeoDict
    from geoist.others.grid2d import Grid2D
    xmin = np.min(x)
    xmax = np.max(x)
    ymin = np.min(y)
    ymax = np.max(y)
    nx = len(x)
    ny = len(y)
    dx = (xmax-xmin)/(nx-1)
    dy = (ymax-ymin)/(ny-1)    
    geodict = {'xmin':xmin,'xmax':xmax,'ymin':ymin,'ymax':ymax,'dx':dx,'dy':dy,'nx':nx,'ny':ny}
    gd = GeoDict(geodict ,adjust = 'res')
    grid = Grid2D(data[::-1], gd)
    return grid
    
def Grid2Xyz(grid):
    """Transform grid data object to x,y,z of 1-D array
    
	Parameters:
		grid is a Grid2D object 
		
    return:
        x,y : The x and y coordinates of the grid points
        z : The value at the grid points
    """   
    xmin,xmax,ymin,ymax = grid.getBounds()
    pdata = grid.getData()[::-1]
    nr,nc = pdata.shape
    xrange = np.linspace(xmin,xmax,num=nc)
    yrange = np.linspace(ymin,ymax,num=nr)
    z = pdata.T
    z = z.ravel()
    x = np.zeros_like(z)
    y = np.zeros_like(z)
    k = 0
    for j in range(len(xrange)):
        for i in range(len(yrange)):
            x[k] = xrange[j] 
            y[k] = yrange[i]
            k += 1
    return x,y,z

def xyz2Grid(x, y ,z):
    """Transform x, y,z to grid data object 
    
	Parameters:
        x,y : The x and y coordinates of the grid points
        z : The value at the grid points     
		  
    return:
        grid :grid is a Grid2D object
    """       
    from geoist.others.geodict import GeoDict
    from geoist.others.grid2d import Grid2D
    xmin = np.min(x)
    xmax = np.max(x)
    ymin = np.min(y)
    ymax = np.max(y)
    nx = len(set(x))
    ny = len(set(y))
    dx = (xmax-xmin)/(nx-1)
    dy = (ymax-ymin)/(ny-1)      
    geodict = {'xmin':xmin,'xmax':xmax,'ymin':ymin,'ymax':ymax,'dx':dx,'dy':dy,'nx':nx,'ny':ny}
    gd = GeoDict(geodict, adjust = 'res')
    data = z.reshape(nx,ny).T[::-1]
    grid = Grid2D(data[::-1], gd)
    return grid

def grid2srf(grid, filename='outsrf.grd', fformat = 'asc'): # bin
    """
     grid is a Grid2D object
    """
    from geoist import DATA_PATH 
    from geoist.pfm.grdio import grddata
    import numpy.ma as ma
    g1out = grddata()
    g1out.cols = grid.getGeoDict().nx
    g1out.rows = grid.getGeoDict().ny
    g1out.xmin = grid.getGeoDict().xmin
    g1out.xmax = grid.getGeoDict().xmax
    g1out.ymin = grid.getGeoDict().ymin
    g1out.ymax = grid.getGeoDict().ymax
    mask  = np.isnan(grid.getData())
    mx = ma.array(grid.getData(),mask=mask)
    g1out.data = mx[::-1] #grid.getData()
    #g1out.data0[np.isnan(grid.getData())] = 1.701410009187828e+38
    if fformat == 'asc':
        g1out.export_surfer(Path(DATA_PATH, filename), True, 'ascii')
    else:
        g1out.export_surfer(Path(DATA_PATH, filename), True, 'binary')
    return g1out
   
def map2DGrid(ax, grid, tstr, xlen=1.0, ylen=1.0, isLeft=False, 
              prj = 'lcc', bous = 0, gmap = 0, xinc = None, yinc = None, 
              pnts = None, cmap='gist_earth'):
    """
    grid is a Grid2D object 
    prj : lcc / merc
    """
    xmin,xmax,ymin,ymax = grid.getBounds()
    pdata = grid.getData()
    nr,nc = pdata.shape
    lonrange = np.linspace(xmin,xmax,num=nc)
    latrange = np.linspace(ymin,ymax,num=nr)
    lon,lat = np.meshgrid(lonrange,latrange)
    
    if xinc != None:
        xmax = xmax + xinc
    if yinc != None:
        ymax = ymax + yinc

    latmean = np.mean([ymin,ymax])
    lonmean = np.mean([xmin,xmax])
    
    #"Lambert Conformal Conic",lcc 'merc': "Mercator",
    if gmap == 1:
        m = Basemap(projection='robin',ellps='WGS84',\
                resolution='c',area_thresh=1000.,\
                lat_0 = latmean, lon_0=lonmean, ax=ax)    
    elif gmap == 2:
        m = Basemap(projection='hammer',ellps='WGS84',\
                resolution='c',area_thresh=1000.,\
                lat_0 = latmean, lon_0=lonmean, ax=ax)            
    else:
        #print(xmin,xmax,ymin,ymax, latmean, lonmean)
        m = Basemap(llcrnrlon=xmin,llcrnrlat=ymin,urcrnrlon=xmax,urcrnrlat=ymax,\
                ellps='WGS84',\
                resolution='c',area_thresh=1000.,projection=prj,\
                lat_0 = latmean, lon_0=lonmean,ax=ax)

        # #print('else')
        # m = Basemap(llcrnrlon=xmin,llcrnrlat=ymin,urcrnrlon=xmax,urcrnrlat=ymax,\
        #         rsphere=(6378137.00,6356752.3142),\
        #         resolution='c',area_thresh=1000.,projection= prj,\
        #         lat_1=latmean,lon_0=lonmean,ax=ax)        
    # draw coastlines and political boundaries.
    m.drawcoastlines()
    if bous == 1:
        m.drawcountries()
    elif bous == 2:
        m.drawcountries()
        m.drawstates()
        
    lons = np.arange(xmin,xmax,xlen)
    lats = np.arange(ymin,ymax,ylen)


    if isLeft:
        labels = labels=[1,0,0,0]
    else:
        labels = labels=[0,0,0,0]
    
    if not gmap:    
        m.drawparallels(lats,labels=labels,color='white',fmt='%.1f') # draw parallels
        m.drawmeridians(lons,labels=[0,0,0,1],color='white',fmt='%.1f') # draw meridians  
    pmesh = m.pcolormesh(lon,lat,np.flipud(grid.getData()),latlon=True, cmap = cmap)
    
    if pnts != None:
        x, y = m(*np.meshgrid(pnts[0], pnts[1])) #(lons,lats))        
        m.scatter(x, y, marker = 'o', color = 'm')

    if ax is not None:
        ax.set_title(tstr)
    m.colorbar(pmesh)
    
def find_nearest(x, x0) -> Tuple[int, Any]:
    """
    This find_nearest function does NOT assume sorted input

    inputs:
    x: array (float, int, datetime, h5py.Dataset) within which to search for x0
    x0: singleton or array of values to search for in x

    outputs:
    idx: index of flattened x nearest to x0  (i.e. works with higher than 1-D arrays also)
    xidx: x[idx]

    Observe how bisect.bisect() gives the incorrect result!

    idea based on:
    http://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array

    """
    x = np.asanyarray(x)  # for indexing upon return
    x0 = np.atleast_1d(x0)

    if x.size == 0 or x0.size == 0:
        raise ValueError('empty input(s)')

    if x0.ndim not in (0, 1):
        raise ValueError('2-D x0 not handled yet')

    ind = np.empty_like(x0, dtype=int)

    # NOTE: not trapping IndexError (all-nan) becaues returning None can surprise with slice indexing
    for i, xi in enumerate(x0):
        if xi is not None and (isinstance(xi, (datetime.datetime, datetime.date, np.datetime64)) or np.isfinite(xi)):
            ind[i] = np.nanargmin(abs(x-xi))
        else:
            raise ValueError('x0 must NOT be None or NaN to avoid surprising None return value')

    return ind.squeeze()[()], x[ind].squeeze()[()]   # [()] to pop scalar from 0d array while being OK with ndim>0

def tickfix(t, fg, ax, tfmt: str='%H:%M:%S'):
    majtick, mintick = timeticks(t[-1] - t[0])
    if majtick:
        ax.xaxis.set_major_locator(majtick)
    if mintick:
        ax.xaxis.set_minor_locator(mintick)
    ax.xaxis.set_major_formatter(DateFormatter(tfmt))
    fg.autofmt_xdate()

    ax.autoscale(True, 'x', tight=True)

    # ax.tick_params(axis='both',which='both')
    # ax.grid(True,which='both')


def timeticks(tdiff):
    """
    NOTE do NOT use "interval" or ticks are misaligned!  use "bysecond" only!
    """
    if isinstance(tdiff, xarray.DataArray):  # len==1
        tdiff = timedelta(seconds=tdiff.values / np.timedelta64(1, 's'))

    assert isinstance(tdiff, timedelta), 'expecting datetime.timedelta'

    if tdiff > timedelta(hours=2):
        return None, None

    elif tdiff > timedelta(minutes=20):
        return MinuteLocator(byminute=range(0, 60, 5)), MinuteLocator(byminute=range(0, 60, 2))

    elif (timedelta(minutes=10) < tdiff) & (tdiff <= timedelta(minutes=20)):
        return MinuteLocator(byminute=range(0, 60, 2)), MinuteLocator(byminute=range(0, 60, 1))

    elif (timedelta(minutes=5) < tdiff) & (tdiff <= timedelta(minutes=10)):
        return MinuteLocator(byminute=range(0, 60, 1)), SecondLocator(bysecond=range(0, 60, 30))

    elif (timedelta(minutes=1) < tdiff) & (tdiff <= timedelta(minutes=5)):
        return SecondLocator(bysecond=range(0, 60, 30)), SecondLocator(bysecond=range(0, 60, 10))

    elif (timedelta(seconds=30) < tdiff) & (tdiff <= timedelta(minutes=1)):
        return SecondLocator(bysecond=range(0, 60, 10)), SecondLocator(bysecond=range(0, 60, 2))

    else:
        return SecondLocator(bysecond=range(0, 60, 2)),  SecondLocator(bysecond=range(0, 60, 1))

# oldfig = plt.figure

# def copyfig(fig=None):
#     # store the image in a buffer using savefig(), this has the
#     # advantage of applying all the default savefig parameters
#     # such as background color; those would be ignored if you simply
#     # grab the canvas
#     if fig is None:
#         # find the figure window that has UI focus right now (not necessarily the same as plt.gcf())
#         fig_window_text = GetWindowText(GetForegroundWindow())
#         for i in plt.get_fignums():
#             if plt.figure(i).canvas.get_window_title() == fig_window_text:
#                 fig = plt.figure(i)
#                 break

#     with io.BytesIO() as buf:
#         fig.savefig(buf)
#         im = Image.open(buf)

#         with io.BytesIO() as output:
#             im.convert("RGB").save(output, "BMP")
#             data = output.getvalue()[14:]  # The file header off-set of BMP is 14 bytes

#     win32clipboard.OpenClipboard()
#     win32clipboard.EmptyClipboard()
#     win32clipboard.SetClipboardData(win32clipboard.CF_DIB, data)  # DIB = device independent bitmap
#     win32clipboard.CloseClipboard()


# def newfig(*args, **kwargs):
#     fig = oldfig(*args, **kwargs)

#     def clipboard_handler(event):
#         if event.key == 'ctrl+c':
#             copyfig()

#     fig.canvas.mpl_connect('key_press_event', clipboard_handler)
#     return fig


# plt.figure = newfig
