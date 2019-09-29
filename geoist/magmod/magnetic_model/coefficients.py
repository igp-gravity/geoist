#-------------------------------------------------------------------------------
#
#  Spherical Harmonic Coefficients.
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
# pylint: disable=too-few-public-methods,abstract-method

from numpy import (
    inf, array, zeros, dot, digitize, argsort, abs as aabs, stack,
)
from ..time_util import decimal_year_to_mjd2000, mjd2000_to_decimal_year


class SHCoefficients(object):
    """ Abstract base class for the spherical harmonic coefficients. """

    def __init__(self, is_internal=True, **kwargs):
        self.is_internal = is_internal
        self.validity = self._get_converted_validity(**kwargs)

    @staticmethod
    def _get_converted_validity(validity_start=None, validity_end=None,
                                to_mjd2000=None, **_):
        """ Get validity range converted to MJD2000 using the optionally
        provided conversion function.
        """

        def _get_converted_value_or_default(value, default, conversion_function):
            if value is None:
                value = default
            if conversion_function is not None:
                value = conversion_function(value)
            return value

        return (
            _get_converted_value_or_default(validity_start, -inf, to_mjd2000),
            _get_converted_value_or_default(validity_end, +inf, to_mjd2000),
        )


    def is_valid(self, time):
        """ Check if the time is within the coefficients validity range. """
        validity_start, validity_end = self.validity
        return (validity_start <= time) & (time <= validity_end)

    @property
    def degree(self):
        """ Get [maximum] model degree. """
        raise NotImplementedError

    @property
    def min_degree(self):
        """ Get minimum model degree.
        Below this degree all model coefficients are zero.
        """
        raise NotImplementedError


    def __call__(self, time, **parameters):
        """ Return the matrix of the full model coefficients. """
        raise NotImplementedError


class CombinedSHCoefficients(SHCoefficients):
    """ Model composed of multiple coefficient sets. """

    def __init__(self, *items):
        if len(items) < 1:
            raise ValueError(
                "The composed model must be composed from at least one "
                "coefficient set."
            )

        item = items[0]
        is_internal = item.is_internal
        validity_start, validity_end = item.validity
        degree = item.degree
        min_degree = item.min_degree

        for item in items[1:]:
            if is_internal != item.is_internal:
                raise ValueError(
                    "Mixing of external and internal coefficient sets!"
                )
            new_start, new_end = item.validity
            validity_start = max(validity_start, new_start)
            validity_end = min(validity_end, new_end)
            degree = max(degree, item.degree)
            min_degree = min(min_degree, item.min_degree)

        SHCoefficients. __init__(
            self, is_internal=is_internal, validity_start=validity_start,
            validity_end=validity_end,
        )

        self._degree = degree
        self._min_degree = min_degree
        self._items = items

    @property
    def degree(self):
        return self._degree

    @property
    def min_degree(self):
        return self._min_degree

    def __call__(self, time, **parameters):
        max_degree = parameters.get("max_degree", -1)
        degree = self.degree if max_degree < 0 else min(self.degree, max_degree)
        coeff_full = zeros((coeff_size(degree), 2))
        for item in self._items:
            item_coeff, item_degree = item(time, **parameters)
            coeff_full[:coeff_size(item_degree), :] += item_coeff
        return coeff_full, degree



class SparseSHCoefficients(SHCoefficients):
    """ Base class for sparse spherical harmonic coefficients. """
    def __init__(self, indices, coefficients, **kwargs):
        SHCoefficients.__init__(self, **kwargs)
        n_idx, m_idx = indices[..., 0], indices[..., 1]
        self._degree = n_idx.max()
        self._min_degree = n_idx.min()
        self._index = stack((
            aabs(m_idx) + (n_idx*(n_idx + 1))//2,
            (m_idx < 0).astype('int'),
            n_idx,
        ), 1)
        self._coeff = coefficients

    @property
    def degree(self):
        return self._degree

    @property
    def min_degree(self):
        return self._min_degree

    def _subset(self, min_degree, max_degree):
        """ Get subset of the coefficients for the give min. and max. degrees.
        """
        default_max_degree = self._degree
        default_min_degree = self._min_degree

        #degree = self._degree
        index = self._index
        coeff = self._coeff

        if max_degree < 0:
            max_degree = default_max_degree

        if min_degree < 0:
            min_degree = default_min_degree

        if min_degree > default_min_degree or max_degree < default_max_degree:
            idx, = (
                (index[:, 2] <= max_degree) & (index[:, 2] >= min_degree)
            ).nonzero()
            coeff = coeff[idx]
            index = index[idx]
            degree = default_max_degree

            if index.shape[0] > 0:
                degree = index[:, 2].max()
            else:
                degree = 0
        else:
            degree = default_max_degree

        return degree, coeff, index[:, 0], index[:, 1]


class SparseSHCoefficientsConstant(SparseSHCoefficients):
    """ Time invariant sparse spherical harmonic coefficients. """
    def __init__(self, indices, coefficients, **kwargs):
        SparseSHCoefficients.__init__(self, indices, coefficients, **kwargs)

    def __call__(self, time, **parameters):
        degree, coeff, index, kind = self._subset(
            parameters.get("min_degree", -1), parameters.get("max_degree", -1)
        )
        coeff_full = zeros((coeff_size(degree), 2))
        coeff_full[index, kind] = coeff
        return coeff_full, degree


class SparseSHCoefficientsTimeDependent(SparseSHCoefficients):
    """ Time dependent sparse spherical harmonic coefficients
    evaluated by piecewise linear interpolation of a time series of
    coefficients snapshots.
    """
    def __init__(self, indices, coefficients, times, **kwargs):
        order = argsort(times)
        kwargs['validity_start'] = kwargs.get('validity_start', times[order[0]])
        kwargs['validity_end'] = kwargs.get('validity_end', times[order[-1]])
        self._times = _convert(times[order], kwargs.get("to_mjd2000"))
        SparseSHCoefficients.__init__(
            self, indices, coefficients[:, order], **kwargs
        )

    def __call__(self, time, **parameters):
        degree, coeff, index, kind = self._subset(
            parameters.get("min_degree", -1), parameters.get("max_degree", -1)
        )
        coeff_full = zeros((coeff_size(degree), 2))
        coeff_full[index, kind] = self._interpolate_coefficients(time, coeff)
        return coeff_full, degree

    def _interpolate_coefficients(self, time, coeff):
        """ Return interpolated coefficients. """
        idx, basis = self._interpolation_basis(time)
        return dot(coeff[:, idx], basis).reshape(coeff.shape[0])

    def _interpolation_basis(self, time):
        """ Return indices and corresponding non-zero basis functions. """
        times = self._times
        idx1 = digitize([time], times)[0]
        idx1 = min(times.size - 1, max(1, idx1))
        idx0 = idx1 - 1
        alpha = (time -times[idx0])/(times[idx1] - times[idx0])
        return [idx0, idx1], array([1.0 - alpha, alpha]).reshape((2, 1))


class SparseSHCoefficientsTimeDependentDecimalYear(SparseSHCoefficientsTimeDependent):
    """ Time dependent sparse spherical harmonic coefficients
    evaluated by piecewise linear interpolation of a time series of
    coefficients snapshots with time in decimal year interpolated
    in the decimal years time domain.
    """

    def __init__(self, indices, coefficients, times,
                 to_mjd2000=decimal_year_to_mjd2000,
                 to_decimal_year=mjd2000_to_decimal_year,
                 **kwargs):
        SparseSHCoefficientsTimeDependent.__init__(
            self, indices, coefficients, times, **kwargs
        )
        self._to_decimal_year = to_decimal_year
        # Fix the validity range to be in the expected MJD2000.
        self.validity = self._get_converted_validity(
            *self.validity, to_mjd2000=to_mjd2000
        )

    def __call__(self, time, **parameters):
        return SparseSHCoefficientsTimeDependent.__call__(
            self, self._to_decimal_year(time), **parameters
        )


def coeff_size(degree):
    """ Size of the full coefficient array. """
    return ((degree + 2)*(degree + 1))//2


def _convert(value, conversion_function):
    """ Convert value using the optional conversion function. """
    if conversion_function is not None:
        value = conversion_function(value)
    return value
