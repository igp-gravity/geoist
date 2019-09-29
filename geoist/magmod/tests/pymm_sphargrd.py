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
from itertools import chain, product
from random import random
from math import pi, cos, sin, sqrt
from numpy import array, zeros
from numpy.testing import assert_allclose
from magmod._pymm import sphargrd
from magmod.tests.pymm_spharpot import SphericalHarmonicsCommonMixIn
from magmod.tests.data import sifm
from magmod.tests.data import mma_external


class SphericalHarmonicsGradientTestMixIn(SphericalHarmonicsCommonMixIn):

    @classmethod
    def reference_gradient(cls, degree, coef_g, coef_h, latitude, longitude, radius):
        rad_series, sin_series, cos_series, p_series, dp_series = cls.get_series(
            degree, latitude, longitude, radius
        )
        m_idx = array(list(chain.from_iterable(
            range(n + 1) for n in range(degree + 1)
        )))
        n_idx = array(list(chain.from_iterable(
            [n]*(n+1)  for n in range(degree + 1)
        )))

        grad_lat = (rad_series[n_idx] * dp_series * (
            coef_g * cos_series[m_idx] + coef_h * sin_series[m_idx]
        )).sum()

        cos_lat = cos(latitude * pi / 180.0)
        if cos_lat > 1e-10:
            grad_lon = -(m_idx * rad_series[n_idx] * p_series * (
                coef_g * sin_series[m_idx] - coef_h * cos_series[m_idx]
            )).sum() / cos_lat
        else:
            sin_lat = sin(latitude * pi / 180.0)
            sin_lon = sin(longitude * pi / 180.0)
            cos_lon = cos(longitude * pi / 180.0)

            scale = [1.0]
            sqn1, ps1, ps0 = 1.0, 1.0, 1.0
            for i in range(2, degree + 1):
                # evaluate ratio between the Gauss-normalised and Schmidt
                # quasi-normalised associated Legendre functions.
                n = float(i)
                tmp = ((n-1)*(n-1)-1)/((2*n-1)*(2*n-3))
                ps1, ps0 = ps0, ps0*sin_lat - ps1*tmp
                sqn1 *= (2*n-1)/n
                scale.append(ps0 * sqn1 * sqrt((n*2)/(n+1)))
            scale = array(scale[:(degree + 1)])
            idx = array([
                1 + (n*(n + 1))//2 for n in range(1, degree + 1)
            ], dtype="int")
            grad_lon = -(scale * rad_series[1:] * (
                coef_g[idx]*sin_lon - coef_h[idx]*cos_lon
            )).sum()

        rad_scale = n_idx + 1 if cls.is_internal else -n_idx
        grad_rad = -(rad_scale * rad_series[n_idx] * p_series * (
            coef_g * cos_series[m_idx] + coef_h * sin_series[m_idx]
        )).sum()

        return array([grad_lat, grad_lon, grad_rad])

    @classmethod
    def eval_gradient(cls, degree, coef_g, coef_h, latitude, longitude, radius):
        rad_series, sin_series, cos_series, p_series, dp_series = cls.get_series(
            degree, latitude, longitude, radius
        )
        return sphargrd(
            latitude, degree, coef_g, coef_h, p_series, dp_series, rad_series,
            sin_series, cos_series, is_internal=cls.is_internal
        )

    def test_coefficients(self):
        max_degree = 3
        coords = [
            (lat, lon, 6371.2) for lat, lon
            in product(range(-90, 91, 10), range(-180, 180, 20))
        ]

        for degree in range(max_degree + 1):
            size = ((degree + 1)*(degree + 2))//2
            offset = (degree*(degree + 1))//2
            for order in range(0, degree + 1):
                coef_g, coef_h = zeros(size), zeros(size)
                coef_g[order + offset] = 1.0

                for latitude, longitude, radius in coords:
                    assert_allclose(
                        self.eval_gradient(
                            degree, coef_g, coef_h, latitude, longitude, radius
                        ),
                        self.reference_gradient(
                            degree, coef_g, coef_h, latitude, longitude, radius
                        ),
                        atol=1e-14
                    )

                coef_g, coef_h = zeros(size), zeros(size)
                coef_h[order + offset] = 1.0

                for latitude, longitude, radius in coords:
                    assert_allclose(
                        self.eval_gradient(
                            degree, coef_g, coef_h, latitude, longitude, radius
                        ),
                        self.reference_gradient(
                            degree, coef_g, coef_h, latitude, longitude, radius
                        ),
                        atol=1e-14
                    )


    def test_gradient(self):
        coords = [
            (lat, lon, 6371.2*(1.0 + random())) for lat, lon
            in product(range(-90, 91, 5), range(-180, 181, 10))
        ]

        for latitude, longitude, radius in coords:
            try:
                assert_allclose(
                    self.eval_gradient(
                        self.degree, self.coef_g, self.coef_h,
                        latitude, longitude, radius
                    ),
                    self.reference_gradient(
                        self.degree, self.coef_g, self.coef_h,
                        latitude, longitude, radius
                    ),
                )
            except AssertionError as exc:
                raise AssertionError(
                    "point coordinates: (%s, %s, %s)\n%s" % (
                        latitude, longitude, radius, str(exc)
                    )
                )

    def test_gradient_compared_with_finite_differences(self):
        # Compare gradient with the calculated finite differences.
        # The test fails if there in no match between the finite difference
        # approximation and the gradient evaluated by the spherical harmonics.
        eps = 0.1
        eps_deg = eps * (180.0 / (pi * 6371.2))
        radius = 6371.2

        def _compare_with_findiff(coord_centre, coord_lat, coord_lon, coord_rad):
            pot0 = self.eval_potential(
                self.degree, self.coef_g, self.coef_h, *coord_centre
            )
            pot_lat = self.eval_potential(
                self.degree, self.coef_g, self.coef_h, *coord_lat
            )
            pot_lon = self.eval_potential(
                self.degree, self.coef_g, self.coef_h, *coord_lon
            )
            pot_rad = self.eval_potential(
                self.degree, self.coef_g, self.coef_h, *coord_rad
            )

            grad_approx = array([
                (pot_lat - pot0)/eps, (pot_lon - pot0)/eps, (pot_rad - pot0)/eps
            ])
            grad_spharm = self.eval_gradient(
                self.degree, self.coef_g, self.coef_h, *coord_centre
            )

            try:
                assert_allclose(grad_spharm, grad_approx, rtol=1e-4, atol=1.0)
            except AssertionError as exc:
                latitude, longitude, radius = coord_centre
                raise AssertionError(
                    "point coordinates: (%s, %s, %s)\n%s" % (
                        latitude, longitude, radius, str(exc)
                    )
                )

        coords = list(product(range(-80, 81, 10), range(-180, 180, 20)))
        for latitude, longitude in coords:
            cos_lat = cos(latitude * pi / 180.0)
            _compare_with_findiff(
                (latitude, longitude, radius),
                (latitude + eps_deg, longitude, radius),
                (latitude, longitude + eps_deg / cos_lat, radius),
                (latitude, longitude, radius + eps)
            )

        for longitude in range(-180, 180, 20):
            _compare_with_findiff(
                (90, longitude, radius),
                (90 - eps_deg, longitude + 180, radius),
                (90 - eps_deg, longitude + 90, radius),
                (90, longitude, radius + eps),
            )

        for longitude in range(-180, 180, 20):
            _compare_with_findiff(
                (-90, longitude, radius),
                (-90 + eps_deg, longitude, radius),
                (-90 + eps_deg, longitude + 90, radius),
                (-90, longitude, radius + eps),
            )


class TestSphericalHarmonicsGradientInternal(TestCase, SphericalHarmonicsGradientTestMixIn):
    is_internal = True
    degree = sifm.DEGREE
    coef_g = sifm.COEF_G
    coef_h = sifm.COEF_H


class TestSphericalHarmonicsGradientExternal(TestCase, SphericalHarmonicsGradientTestMixIn):
    is_internal = False
    degree = mma_external.DEGREE
    coef_g = mma_external.COEF_Q
    coef_h = mma_external.COEF_S


if __name__ == "__main__":
    main()
