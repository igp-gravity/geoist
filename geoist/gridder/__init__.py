"""
Create and operate on data grids, scatters, and profiles.
"""

from .slicing import inside, cut
from .interpolation import interp, interp_at, profile
from .padding import pad_array, unpad_array, pad_coords
from .genpnt import regular, scatter, circular_scatter

