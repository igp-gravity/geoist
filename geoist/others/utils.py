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
