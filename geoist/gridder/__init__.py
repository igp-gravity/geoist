"""
Create and operate on data grids, scatters, and profiles.
"""
import glob
from os.path import dirname,basename,isfile
modules = glob.glob(dirname(__file__)+"*.py")
__all__ = [basename(f)[:-3] for f in modules if isfile(f) and not f.endswith('__init__.py')]

from .slicing import inside, cut
from .interpolation import interp, interp_at, profile
from .padding import pad_array, unpad_array, pad_coords
from .point_generation import regular, scatter, circular_scatter
from .utils import spacing
