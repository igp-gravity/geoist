"""
The ``fatiando`` package contains subpackages for different geophysical methods
and other utilities.

See the API reference for each subpackage for a list of all functions and
classes defined by it.
"""

from ._version import get_versions
import glob
from os.path import dirname,basename,isfile
modules = glob.glob(dirname(__file__)+"*.py")
__all__ = [basename(f)[:-3] for f in modules if isfile(f) and not f.endswith('__init__.py')]


__version__ = get_versions()['version']
__commit__ = get_versions()['full']
del get_versions


def test(doctest=True, verbose=False, coverage=False):
    """
    Run the test suite for Fatiando a Terra.

    Uses `py.test <http://pytest.org/>`__ to discover and run the tests. If you
    haven't already, you can install it with `conda
    <http://conda.pydata.org/>`__ or `pip <https://pip.pypa.io/en/stable/>`__.

    Parameters:

    * doctest : bool
        If ``True``, will run the doctests as well (code examples that start
        with a ``>>>`` in the docs).
    * verbose : bool
        If ``True``, will print extra information during the test run.
    * coverage : bool
        If ``True``, will run test coverage analysis on the code as well.
        Requires ``pytest-cov``.

    Raises:

    * ``AssertionError`` if pytest returns a non-zero error code indicating
      that some tests have failed.

    """
    import pytest
    args = []
    if verbose:
        args.append('-v')
    if coverage:
        args.append('--cov=geoist')
        args.append('--cov-report=term-missing')
    if doctest:
        args.append('--doctest-modules')
    args.append('--pyargs')
    args.append('geoist')
    status = pytest.main(args)
    assert status == 0, "Some tests have failed."

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
