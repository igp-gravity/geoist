#-------------------------------------------------------------------------------
#
#  Swarm MIO_SHA_2* product spherical harmonic coefficients.
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
from io import open
from numpy import inf, nan
from numpy.testing import assert_allclose
from geoist.magmod.time_util import decimal_year_to_mjd2000
from geoist.magmod.magnetic_time import mjd2000_to_magnetic_universal_time
from geoist.magmod.magnetic_model.coefficients_mio import SparseSHCoefficientsMIO
from geoist.magmod.magnetic_model.tests.data import SWARM_MIO_SHA_2_TEST_DATA
from geoist.magmod.magnetic_model.parser_mio import parse_swarm_mio_file


class MIOSHCoeffMixIn(object):
    is_internal = None
    degree = 2
    min_degree = 1
    validity = (-inf, inf)

    def test_min_degree(self):
        self.assertEqual(self.coefficients.min_degree, self.min_degree)

    def test_degree(self):
        self.assertEqual(self.coefficients.degree, self.degree)

    def test_is_internal(self):
        self.assertEqual(self.coefficients.is_internal, self.is_internal)

    def test_validity(self):
        assert_allclose(self.coefficients.validity, self.validity)

    def test_is_valid_success(self):
        self.assertTrue(self.coefficients.is_valid(0.0))

    def test_is_valid_fail_nan(self):
        self.assertFalse(self.coefficients.is_valid(nan))

    def eval_coeff(self, time, lat_sol=None, lon_sol=None, **options):
        return self.coefficients(time, mjd2000_to_magnetic_universal_time(
            time, lat_ngp=self.lat_ngp, lon_ngp=self.lon_ngp,
            lat_sol=lat_sol, lon_sol=lon_sol,
        ), **options)


class TestSparseSHCoefficientsMIOInternal(TestCase, MIOSHCoeffMixIn):
    is_internal = True

    def setUp(self):
        with open(SWARM_MIO_SHA_2_TEST_DATA, encoding="ascii") as file_in:
            data = parse_swarm_mio_file(file_in)

        self.lat_ngp = data["lat_NGP"]
        self.lon_ngp = data["lon_NGP"]
        self.coefficients = SparseSHCoefficientsMIO(
            data["nm"], data["gh"],
            ps_extent=(data["pmin"], data["pmax"], data["smin"], data["smax"]),
            is_internal=True,
        )

    def test_callable(self):
        coeff, degree = self.eval_coeff(decimal_year_to_mjd2000(2018.5))
        assert_allclose(coeff, [
            (0, 0),
            (0.11415226, 0), (-0.22361224, -0.02624283),
            (-0.23244221, 0), (-0.34379927, -0.81524805), (0.01412196, -0.18049787),
        ], atol=1e-8)
        self.assertEqual(degree, self.degree)

    def test_extra_sub_solar_point(self):
        coeff, degree = self.eval_coeff(
            decimal_year_to_mjd2000(2018.5), lat_sol=0.0, lon_sol=0.0,
        )
        assert_allclose(coeff, [
            (0.0, 0.0),
            (0.12120376078924046, 0.0),
            (-0.22660807757361834, 0.008720797143766132),
            (-0.2485487953382298, 0.0),
            (-0.40755502000159416, -0.7608326017305439),
            (-0.04181600178841578, -0.16655581670562275),
        ], atol=1e-8)
        self.assertEqual(degree, self.degree)

    def test_callable_max_degree(self):
        coeff, degree = self.eval_coeff(
            decimal_year_to_mjd2000(2018.5), max_degree=1
        )
        assert_allclose(coeff, [
            (0, 0),
            (0.11415226, 0), (-0.22361224, -0.02624283),
        ], atol=1e-8)
        self.assertEqual(degree, 1)

    def test_callable_min_degree(self):
        coeff, degree = self.eval_coeff(
            decimal_year_to_mjd2000(2018.5), min_degree=2
        )
        assert_allclose(coeff, [
            (0, 0),
            (0, 0), (0, 0),
            (-0.23244221, 0), (-0.34379927, -0.81524805), (0.01412196, -0.18049787),
        ], atol=1e-8)
        self.assertEqual(degree, self.degree)

    def test_callable_min_max_degree(self):
        coeff, degree = self.eval_coeff(
            decimal_year_to_mjd2000(2018.5), max_degree=1, min_degree=1
        )
        assert_allclose(coeff, [
            (0, 0),
            (0.11415226, 0), (-0.22361224, -0.02624283),
        ], atol=1e-8)
        self.assertEqual(degree, 1)

    def test_callable_all_zero(self):
        coeff, degree = self.eval_coeff(
            decimal_year_to_mjd2000(2018.5), min_degree=3
        )
        assert_allclose(coeff, [
            [0, 0],
        ])
        self.assertEqual(degree, 0)


class TestSparseSHCoefficientsMIOExternal(TestCase, MIOSHCoeffMixIn):
    is_internal = False

    def setUp(self):
        with open(SWARM_MIO_SHA_2_TEST_DATA, encoding="ascii") as file_in:
            data = parse_swarm_mio_file(file_in)

        self.lat_ngp = data["lat_NGP"]
        self.lon_ngp = data["lon_NGP"]
        self.coefficients = SparseSHCoefficientsMIO(
            data["nm"], data["qs"],
            ps_extent=(data["pmin"], data["pmax"], data["smin"], data["smax"]),
            is_internal=False,
        )

    def test_callable(self):
        coeff, degree = self.eval_coeff(decimal_year_to_mjd2000(2018.5))
        assert_allclose(coeff, [
            (0, 0),
            (-1.44080531, 0), (-0.27711333, -0.25918073),
            (-0.63645813, 0), (-0.22583255, -2.02305245), (0.05941826, -0.24358544),
        ], atol=1e-8)
        self.assertEqual(degree, self.degree)


if __name__ == "__main__":
    main()
