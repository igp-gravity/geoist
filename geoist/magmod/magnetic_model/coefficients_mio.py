#-------------------------------------------------------------------------------
#
#  Spherical Harmonic Coefficients specific to Swarm MIO model.
#
# Author: Martin Paces <martin.paces@eox.at>
#
#-------------------------------------------------------------------------------
# Copyright (C) 2018 EOX IT Services GmbH
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies of this Software or works derived from this Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#-------------------------------------------------------------------------------
#pylint: disable=no-name-in-module

from math import pi
from numpy import zeros, arange, broadcast_to, sin, cos
from .coefficients import SparseSHCoefficients, coeff_size
from ..magnetic_time import mjd2000_to_magnetic_universal_time
from ..time_util import mjd2000_to_year_fraction

F_SEASONAL = 2*pi
F_DIURNAL = 2*pi/24.


class SparseSHCoefficientsMIO(SparseSHCoefficients):
    """ Time dependent 2D Fourier series Swarm MIO coefficients.
        Parameters:
            indices - array if the nm indices
            coefficients - gh or qs coefficients
            ps_extent - (pmin, pmax, smin, smax) tuple
            is_internal - set False for an external model
    """

    def __init__(self, indices, coefficients, ps_extent,
                 mjd2000_to_year_fraction=mjd2000_to_year_fraction, **kwargs):
        SparseSHCoefficients.__init__(self, indices, coefficients, **kwargs)
        pmin, pmax, smin, smax = ps_extent
        if pmin > pmax or smin > smax:
            raise Exception("Invalid ps_extent %s!" % ps_extent)
        self.ps_extent = (pmin, pmax, smin, smax)
        self.mjd2000_to_year_fraction = mjd2000_to_year_fraction

    def __call__(self, time, mut, **parameters):
        degree, coeff, index, kind = self._subset(
            parameters.get("min_degree", -1), parameters.get("max_degree", -1)
        )
        coeff_full = zeros((coeff_size(degree), 2))
        coeff_full[index, kind] = self._eval_coeff_fourier2d(coeff, time, mut)
        return coeff_full, degree

    def _eval_coeff_fourier2d(self, coeff, mjd2000, mut):
        """ Evaluate model coefficients using the 2D Fourier series. """
        sin_f, cos_f = self._get_sincos_matrices(coeff.shape[0], mjd2000, mut)
        return (coeff[..., 0]*cos_f + coeff[..., 1]*sin_f).sum(axis=(1, 2))

    def _get_sincos_matrices(self, n_coeff, mjd2000, mut):
        """ Get sin/cos matrices used by the 2D Fourier transform. """
        pmin, pmax, smin, smax = self.ps_extent
        n_col = pmax - pmin + 1
        n_row = smax - smin + 1
        f0_seasonal = F_SEASONAL * self.mjd2000_to_year_fraction(mjd2000)
        f0_diurnal = F_DIURNAL * mut
        f_diurnal = f0_diurnal * arange(pmin, pmax + 1)
        f_seasonal = f0_seasonal * arange(smin, smax + 1)
        f_combined = (
            broadcast_to(f_diurnal, (n_row, n_col)) +
            broadcast_to(f_seasonal, (n_col, n_row)).transpose()
        )
        sin_f = broadcast_to(sin(f_combined), (n_coeff, n_row, n_col))
        cos_f = broadcast_to(cos(f_combined), (n_coeff, n_row, n_col))
        return sin_f, cos_f
