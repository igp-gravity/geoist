#-------------------------------------------------------------------------------
#
#  Spherical Harmonic Expansion - Geomagnetic Model - tests
#
#
# Author: Steve Shi Chen <chenshi80@gmail.com>
# 
# Original Author: Martin Paces <martin.paces@eox.at>
#-------------------------------------------------------------------------------
# Copyright (C) 2019 Geoist team
#
#-------------------------------------------------------------------------------
# pylint: disable=missing-docstring

from unittest import TestCase, main
from math import pi, sin, sqrt
from numpy import array
from numpy.testing import assert_allclose
from scipy.special import legendre as legendre_polynomial
from magmod._pymm import legendre


class TestLegendreFunctions(TestCase):

    @staticmethod
    def _assert_allclose(result0, result1):
        p_series0, dp_series0 = result0
        p_series1, dp_series1 = result1
        assert_allclose(dp_series0, dp_series1, atol=1e-6)
        assert_allclose(p_series0, p_series1, atol=1e-6)

    @classmethod
    def reference(cls, latitude, degree):
        return cls.schmidt_norm_legendre(sin(latitude * pi / 180.0), degree)

    @staticmethod
    def schmidt_norm_legendre(x, degree):
        """ Evaluate Schmidt semi-normalised Legendre functions and their
        derivatives.
        """
        # pylint: disable=invalid-name

        # TODO: Get more robust algorithm.
        # NOTE: This simple reference implementation is prone to truncation
        # errors and overflows for higher degrees.

        def _eval_functions():
            y = (1.0 - x*x)
            yield legendre_polynomial(0)(x)
            for n in range(1, degree + 1):
                leg_pol = legendre_polynomial(n)
                yield leg_pol(x)
                scale = sqrt(y / ((n * (n + 1))//2))
                yield scale * leg_pol.deriv(1)(x)
                for m  in range(2, n + 1):
                    scale *= sqrt(y / ((n + m) * (n - m + 1)))
                    yield scale * leg_pol.deriv(m)(x)

        def _eval_derivatives():
            yield 0.0
            yield p_series[2]
            yield -p_series[1]
            for n in range(2, degree + 1):
                offset = (n*(n + 1))//2
                yield sqrt((n*(n + 1))//2) * p_series[offset + 1]
                yield 0.5 * (
                    sqrt((n+2)*(n-1)) * p_series[offset + 2] -
                    sqrt(2*n*(n + 1)) * p_series[offset]
                )
                for m in range(2, n):
                    yield 0.5 * (
                        sqrt((n+m+1)*(n-m)) * p_series[offset + m + 1] -
                        sqrt((n+m)*(n-m+1)) * p_series[offset + m - 1]
                    )
                yield -sqrt(n/2.0) * p_series[offset + n - 1]

        p_series = list(_eval_functions())
        dp_series = list(_eval_derivatives())

        return array(p_series), array(dp_series)

    def test_legendre_functions(self):
        degree = 20
        latitudes = [float(v) for v in range(-90, 91, 5)]
        for latitude in latitudes:
            self._assert_allclose(
                legendre(latitude, degree), self.reference(latitude, degree)
            )

    def test_legendre_functions_zero_degree(self):
        self._assert_allclose(legendre(0, 0), ([1.0], [0.0]))
        self._assert_allclose(legendre(-90, 0), ([1.0], [0.0]))
        self._assert_allclose(legendre(+90, 0), ([1.0], [0.0]))

    def test_invalid_degree(self):
        self.assertRaises(ValueError, legendre, 0.0, -1)

    def test_invalid_latitude(self):
        self.assertRaises(ValueError, legendre, -91.0, 1)
        self.assertRaises(ValueError, legendre, +91.0, 1)


if __name__ == "__main__":
    main()
