#-------------------------------------------------------------------------------
#
#  Spherical Harmonic Expansion - Geomagnetic Model - tests
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
from math import pi, sin, cos
from numpy import array
from numpy.testing import assert_allclose
from magmod._pymm import lonsincos


class TestLongitudialSinCosSeries(TestCase):

    @staticmethod
    def reference(value, degree):
        """ Reference implementation. """
        value *= pi / 180.0
        return (
            array([sin(i*value) for i in range(degree + 1)]),
            array([cos(i*value) for i in range(degree + 1)])
        )

    @staticmethod
    def _assert_allclose(result0, result1):
        sin_series0, cos_series0 = result0
        sin_series1, cos_series1 = result1
        assert_allclose(cos_series0, cos_series1, atol=1e-14)
        assert_allclose(sin_series0, sin_series1, atol=1e-14)

    def _test_lonsincos(self, *args, **kwargs):
        degree = 6
        longitudes = [float(v) for v in range(-180, 180, 30)]

        for longitude in longitudes:
            self._assert_allclose(
                lonsincos(longitude, degree, *args, **kwargs),
                self.reference(longitude, degree)
            )

    def test_invalid_degree(self):
        self.assertRaises(ValueError, lonsincos, 0.0, -1)

    def test_lonsincos_zero_degree(self):
        self._assert_allclose(lonsincos(0, 0), ([0.0], [1.0]))
        self._assert_allclose(lonsincos(90, 0, False), ([0.0], [1.0]))
        self._assert_allclose(lonsincos(-90, 0, True), ([0.0], [1.0]))

    def test_lonsincos_default(self):
        self._test_lonsincos()

    def test_lonsincos_fast(self):
        self._test_lonsincos(True)

    def test_lonsincos_slow(self):
        self._test_lonsincos(False)


if __name__ == "__main__":
    main()
