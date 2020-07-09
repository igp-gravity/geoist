#-------------------------------------------------------------------------------
#
#  Spherical Harmonic Coefficients.
#
# Author: Steve Shi Chen <chenshi80@gmail.com>
# 
# Original Author: Martin Paces <martin.paces@eox.at>
#-------------------------------------------------------------------------------
# Copyright (C) 2019 Geoist team
#
#-------------------------------------------------------------------------------
# pylint: disable=missing-docstring, line-too-long, too-few-public-methods

from unittest import TestCase, main
from numpy import nan, inf, isinf, array, dot
from numpy.testing import assert_allclose
from geoist.magmod import decimal_year_to_mjd2000
from geoist.magmod.magnetic_model.coefficients import (
    SparseSHCoefficientsTimeDependent,
    SparseSHCoefficientsTimeDependentDecimalYear,
    SparseSHCoefficientsConstant,
    CombinedSHCoefficients,
)


class SHCoefficinetTestMixIn(object):

    def test_min_degree(self):
        self.assertEqual(self.coefficients.min_degree, self.min_degree)

    def test_degree(self):
        self.assertEqual(self.coefficients.degree, self.degree)

    def test_is_internal(self):
        self.assertEqual(self.coefficients.is_internal, self.is_internal)

    def test_validity(self):
        self.assertEqual(self.coefficients.validity, self.validity)

    def test_is_valid_success(self):
        validity_start, validity_end = self.coefficients.validity
        if isinf(validity_start):
            if isinf(validity_end):
                time = 0.0
            else:
                time = validity_end - 1.0
        else:
            if isinf(validity_end):
                time = validity_start + 1.0
            else:
                time = 0.5*(validity_start + validity_end)
        self.assertTrue(self.coefficients.is_valid(time))

    def test_is_valid_fail_before(self):
        validity_start, _ = self.coefficients.validity
        if not isinf(validity_start):
            self.assertFalse(self.coefficients.is_valid(validity_start - 1.0))

    def test_is_valid_fail_after(self):
        _, validity_end = self.coefficients.validity
        if not isinf(validity_end):
            self.assertFalse(self.coefficients.is_valid(validity_end + 1.0))

    def test_is_valid_fail_nan(self):
        self.assertFalse(self.coefficients.is_valid(nan))

#-------------------------------------------------------------------------------

class CombinedSHCoefficientsMixIn(object):
    is_internal = True
    times0 = array([2012.0, 2016.0, 2014.0])
    indices0 = array([(1, 0), (1, 1), (1, -1)])
    coeff0 = array([
        [1, 3, 2],
        [5, 15, 10],
        [10, 30, 20],
    ])
    indices1 = array([(2, 0), (2, 1), (2, -1), (2, 2), (2, -2)])
    coeff1 = array([1, 5, 10, 8, 12])
    options0 = {}
    options1 = {}
    degree = 2
    min_degree = 1
    validity = (2012.0, 2016.0)

    @property
    def coefficients(self):
        return CombinedSHCoefficients(
            SparseSHCoefficientsTimeDependent(
                self.indices0, self.coeff0, self.times0, **self.options0
            ),
            SparseSHCoefficientsConstant(
                self.indices1, self.coeff1, **self.options1
            )
        )


class TestCombinedSHCoefficientsDefault(TestCase, SHCoefficinetTestMixIn, CombinedSHCoefficientsMixIn):
    def test_callable(self):
        coeff, degree = self.coefficients(2013.0)
        assert_allclose(coeff, [
            [0., 0.], [1.5, 0], [7.5, 15.0], [1, 0], [5, 10.0], [8., 12.],
        ])
        self.assertEqual(degree, self.degree)


class TestCombinedSHCoefficientsMinDegree(TestCase, SHCoefficinetTestMixIn, CombinedSHCoefficientsMixIn):
    def test_callable(self):
        coeff, degree = self.coefficients(2013.0, min_degree=2)
        assert_allclose(coeff, [
            [0., 0.], [0., 0.], [0., 0.], [1., 0.], [5, 10.0], [8., 12.],
        ])
        self.assertEqual(degree, self.degree)


class TestCombinedSHCoefficientsMaxDegree(TestCase, SHCoefficinetTestMixIn, CombinedSHCoefficientsMixIn):
    def test_callable(self):
        coeff, degree = self.coefficients(2013.0, max_degree=1)
        assert_allclose(coeff, [[0., 0.], [1.5, 0], [7.5, 15.0]])
        self.assertEqual(degree, 1)


class TestCombinedSHCoefficientsMinMaxDegree(TestCase, SHCoefficinetTestMixIn, CombinedSHCoefficientsMixIn):
    def test_callable(self):
        coeff, degree = self.coefficients(2013.0, min_degree=1, max_degree=1)
        assert_allclose(coeff, [[0., 0.], [1.5, 0], [7.5, 15.0]])
        self.assertEqual(degree, 1)


class TestCombinedSHCoefficientsInternal(TestCombinedSHCoefficientsDefault):
    is_internal = True
    options0 = {"is_internal": True}
    options1 = {"is_internal": True}


class TestCombinedSHCoefficientsExternal(TestCombinedSHCoefficientsDefault):
    is_internal = False
    options0 = {"is_internal": False}
    options1 = {"is_internal": False}


class TestCombinedSHCoefficientsMixed(TestCase, CombinedSHCoefficientsMixIn):
    options0 = {"is_internal": True}
    options1 = {"is_internal": False}

    def test_mixed_type_failure(self):
        with self.assertRaises(ValueError):
            self.coefficients # pylint: disable=pointless-statement

#-------------------------------------------------------------------------------

class TestSparseSHCoefficientsConstantDefault(TestCase, SHCoefficinetTestMixIn):
    indices = array([(1, 0), (1, 1), (1, -1)])
    coeff = array([1, 5, 10])
    degree = 1
    min_degree = 1
    is_internal = True
    validity = (-inf, inf)
    options = {}

    @property
    def coefficients(self):
        return SparseSHCoefficientsConstant(
            self.indices, self.coeff, **self.options
        )

    def test_callable(self):
        coeff, degree = self.coefficients(2013.0)
        assert_allclose(coeff, [[0., 0.], [1, 0], [5, 10.0]])
        self.assertEqual(degree, self.degree)


class TestSparseSHCoefficientsConstantInternal(TestSparseSHCoefficientsConstantDefault):
    is_internal = True
    options = {"is_internal": True}


class TestSparseSHCoefficientsConstantExternal(TestSparseSHCoefficientsConstantDefault):
    is_internal = False
    options = {"is_internal": False}


class TestSparseSHCoefficientsConstantExtendedValidity(TestSparseSHCoefficientsConstantDefault):
    validity = (2010.0, 2018.0)
    options = {"validity_start": 2010.0, "validity_end": 2018.0}

#-------------------------------------------------------------------------------

class TestSparseSHCoefficientsTimeDependentDefault(TestCase, SHCoefficinetTestMixIn):
    times = array([2012.0, 2016.0, 2014.0])
    indices = array([(1, 0), (1, 1), (1, -1)])
    coeff = array([
        [1, 3, 2],
        [5, 15, 10],
        [10, 30, 20],
    ])
    degree = 1
    min_degree = 1
    is_internal = True
    validity = (2012.0, 2016.0)
    options = {}

    @property
    def coefficients(self):
        return SparseSHCoefficientsTimeDependent(
            self.indices, self.coeff, self.times, **self.options
        )

    def test_callable(self):
        coeff, degree = self.coefficients(2013.0)
        assert_allclose(coeff, [[0., 0.], [1.5, 0], [7.5, 15.0]])
        self.assertEqual(degree, self.degree)

    def test_callable_before_first_time(self):
        coeff, degree = self.coefficients(2011.0)
        assert_allclose(coeff, [[0., 0.], [0.5, 0], [2.5, 5.0]])
        self.assertEqual(degree, self.degree)


class TestSparseSHCoefficientsTimeDependentInternal(TestSparseSHCoefficientsTimeDependentDefault):
    is_internal = True
    options = {"is_internal": True}


class TestSparseSHCoefficientsTimeDependentExternal(TestSparseSHCoefficientsTimeDependentDefault):
    is_internal = False
    options = {"is_internal": False}


class TestSparseSHCoefficientsTimeDependentExtendedValidity(TestSparseSHCoefficientsTimeDependentDefault):
    validity = (2010.0, 2018.0)
    options = {"validity_start": 2010.0, "validity_end": 2018.0}

#-------------------------------------------------------------------------------


class TestSparseSHCoefficientsTimeDependentConvertedDefault(TestCase, SHCoefficinetTestMixIn):
    times = array([2012.0, 2016.0, 2014.0])
    indices = array([(1, 0), (1, 1), (1, -1)])
    coeff = array([
        [1, 3, 2],
        [5, 15, 10],
        [10, 30, 20],
    ])
    degree = 1
    min_degree = 1
    is_internal = True
    validity = tuple(decimal_year_to_mjd2000((2012.0, 2016.0)))
    options = {}

    @property
    def coefficients(self):
        return SparseSHCoefficientsTimeDependent(
            self.indices, self.coeff, self.times,
            to_mjd2000=decimal_year_to_mjd2000, **self.options
        )

    def test_callable(self):
        coeff, degree = self.coefficients(
            dot((0.5, 0.5), decimal_year_to_mjd2000((2012.0, 2014.0)))
        )
        assert_allclose(coeff, [[0., 0.], [1.5, 0], [7.5, 15.0]])
        self.assertEqual(degree, self.degree)

    def test_callable_before_first_time(self):
        coeff, degree = self.coefficients(
            dot((1.5, -0.5), decimal_year_to_mjd2000((2012.0, 2014.0)))
        )
        assert_allclose(coeff, [[0., 0.], [0.5, 0], [2.5, 5.0]])
        self.assertEqual(degree, self.degree)


class TestSparseSHCoefficientsTimeDependentConvertedInternal(TestSparseSHCoefficientsTimeDependentConvertedDefault):
    is_internal = True
    options = {"is_internal": True}


class TestSparseSHCoefficientsTimeDependentConvertedExternal(TestSparseSHCoefficientsTimeDependentConvertedDefault):
    is_internal = False
    options = {"is_internal": False}


class TestSparseSHCoefficientsTimeDependentConvertedExtendedValidity(TestSparseSHCoefficientsTimeDependentConvertedDefault):
    validity = tuple(decimal_year_to_mjd2000((2010.0, 2018.0)))
    options = {"validity_start": 2010.0, "validity_end": 2018.0}

#-------------------------------------------------------------------------------

class TestSparseSHCoefficientsTimeDependentDecimalYearDefault(TestCase, SHCoefficinetTestMixIn):
    times = array([2012.0, 2016.0, 2014.0])
    indices = array([(1, 0), (1, 1), (1, -1)])
    coeff = array([
        [1, 3, 2],
        [5, 15, 10],
        [10, 30, 20],
    ])
    degree = 1
    min_degree = 1
    is_internal = True
    validity = tuple(decimal_year_to_mjd2000((2012.0, 2016.0)))
    options = {}

    @property
    def coefficients(self):
        return SparseSHCoefficientsTimeDependentDecimalYear(
            self.indices, self.coeff, self.times, **self.options
        )

    def test_callable(self):
        coeff, degree = self.coefficients(decimal_year_to_mjd2000(2013.0))
        assert_allclose(coeff, [[0., 0.], [1.5, 0], [7.5, 15.0]])
        self.assertEqual(degree, self.degree)

    def test_callable_before_first_time(self):
        coeff, degree = self.coefficients(decimal_year_to_mjd2000(2011.0))
        assert_allclose(coeff, [[0., 0.], [0.5, 0], [2.5, 5.0]])
        self.assertEqual(degree, self.degree)


class TestSparseSHCoefficientsTimeDependentDecimalYearInternal(TestSparseSHCoefficientsTimeDependentDecimalYearDefault):
    is_internal = True
    options = {"is_internal": True}


class TestSparseSHCoefficientsTimeDependentDecimalYearExternal(TestSparseSHCoefficientsTimeDependentDecimalYearDefault):
    is_internal = False
    options = {"is_internal": False}


class TestSparseSHCoefficientsTimeDependentDecimalYearExtendedValidity(TestSparseSHCoefficientsTimeDependentDecimalYearDefault):
    validity = tuple(decimal_year_to_mjd2000((2010.0, 2018.0)))
    options = {"validity_start": 2010.0, "validity_end": 2018.0}

#-------------------------------------------------------------------------------

class TestSparseSHCoefficientsConstantSubset(TestCase, SHCoefficinetTestMixIn):
    indices = array([(1, 0), (2, 1), (3, -2)])
    coeff = array([1, 5, 10])
    degree = 3
    min_degree = 1
    is_internal = True
    validity = (-inf, inf)
    options = {}

    @property
    def coefficients(self):
        return SparseSHCoefficientsConstant(
            self.indices, self.coeff, **self.options
        )

    def test_callable(self):
        coeff, degree = self.coefficients(2013.0)
        assert_allclose(coeff, [
            [0, 0],
            [1, 0], [0, 0],
            [0, 0], [5, 0], [0, 0],
            [0, 0], [0, 0], [0, 10], [0, 0],
        ])
        self.assertEqual(degree, 3)

    def test_callable_min_degree(self):
        coeff, degree = self.coefficients(2013.0, min_degree=2)
        assert_allclose(coeff, [
            [0, 0],
            [0, 0], [0, 0],
            [0, 0], [5, 0], [0, 0],
            [0, 0], [0, 0], [0, 10], [0, 0],
        ])
        self.assertEqual(degree, 3)

    def test_callable_max_degree(self):
        coeff, degree = self.coefficients(2013.0, max_degree=2)
        assert_allclose(coeff, [
            [0, 0],
            [1, 0], [0, 0],
            [0, 0], [5, 0], [0, 0],
        ])
        self.assertEqual(degree, 2)

    def test_callable_min_max_degree(self):
        coeff, degree = self.coefficients(2013.0, min_degree=2, max_degree=2)
        assert_allclose(coeff, [
            [0, 0],
            [0, 0], [0, 0],
            [0, 0], [5, 0], [0, 0],
        ])
        self.assertEqual(degree, 2)

    def test_callable_all_zero(self):
        coeff, degree = self.coefficients(2013.0, min_degree=4)
        assert_allclose(coeff, [
            [0, 0],
        ])
        self.assertEqual(degree, 0)


class TestSparseSHCoefficientsTimeDependentSubset(TestCase, SHCoefficinetTestMixIn):
    times = array([2012.0, 2014.0])
    indices = array([(1, 0), (2, 1), (3, -2)])
    coeff = array([
        [1, 2],
        [5, 10],
        [10, 20],
    ])
    degree = 3
    min_degree = 1
    is_internal = True
    validity = (2012.0, 2014.0)
    options = {}

    @property
    def coefficients(self):
        return SparseSHCoefficientsTimeDependent(
            self.indices, self.coeff, self.times, **self.options
        )

    def test_callable(self):
        coeff, degree = self.coefficients(2013.0)
        assert_allclose(coeff, [
            [0, 0],
            [1.5, 0], [0, 0],
            [0, 0], [7.5, 0], [0, 0],
            [0, 0], [0, 0], [0, 15], [0, 0],
        ])
        self.assertEqual(degree, 3)

    def test_callable_min_degree(self):
        coeff, degree = self.coefficients(2013.0, min_degree=2)
        assert_allclose(coeff, [
            [0, 0],
            [0, 0], [0, 0],
            [0, 0], [7.5, 0], [0, 0],
            [0, 0], [0, 0], [0, 15], [0, 0],
        ])
        self.assertEqual(degree, 3)

    def test_callable_max_degree(self):
        coeff, degree = self.coefficients(2013.0, max_degree=2)
        assert_allclose(coeff, [
            [0, 0],
            [1.5, 0], [0, 0],
            [0, 0], [7.5, 0], [0, 0],
        ])
        self.assertEqual(degree, 2)

    def test_callable_min_max_degree(self):
        coeff, degree = self.coefficients(2013.0, min_degree=2, max_degree=2)
        assert_allclose(coeff, [
            [0, 0],
            [0, 0], [0, 0],
            [0, 0], [7.5, 0], [0, 0],
        ])
        self.assertEqual(degree, 2)

    def test_callable_all_zero(self):
        coeff, degree = self.coefficients(2013.0, min_degree=4)
        assert_allclose(coeff, [
            [0, 0],
        ])
        self.assertEqual(degree, 0)

#-------------------------------------------------------------------------------

if __name__ == "__main__":
    main()
