from typing import Tuple, Any
import numpy as np
import datetime


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
# %%
    if x.size == 0 or x0.size == 0:
        raise ValueError('empty input(s)')

    if x0.ndim not in (0, 1):
        raise ValueError('2-D x0 not handled yet')
# %%
    ind = np.empty_like(x0, dtype=int)

    # NOTE: not trapping IndexError (all-nan) becaues returning None can surprise with slice indexing
    for i, xi in enumerate(x0):
        if xi is not None and (isinstance(xi, (datetime.datetime, datetime.date, np.datetime64)) or np.isfinite(xi)):
            ind[i] = np.nanargmin(abs(x-xi))
        else:
            raise ValueError('x0 must NOT be None or NaN to avoid surprising None return value')

    return ind.squeeze()[()], x[ind].squeeze()[()]   # [()] to pop scalar from 0d array while being OK with ndim>0
