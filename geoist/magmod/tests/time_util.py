#-------------------------------------------------------------------------------
#
#  Time conversion utilities - test
#
# Author: Steve Shi Chen <chenshi80@gmail.com>
# 
# Original Author: Martin Paces <martin.paces@eox.at>
#-------------------------------------------------------------------------------
# Copyright (C) 2019 Geoist team
#
#-------------------------------------------------------------------------------
# pylint: disable=missing-docstring, invalid-name, too-few-public-methods

from unittest import TestCase, main
from numpy import vectorize, inf, nan
from numpy.random import uniform
from numpy.testing import assert_allclose
from geoist.magmod.time_util import (
    decimal_year_to_mjd2000_simple,
    mjd2000_to_decimal_year_simple,
    mjd2000_to_year_fraction_simple,
)


class TestMjd2000ToYearFractionSimple(TestCase):

    @staticmethod
    def reference(value):
        return vectorize(_mjd2000_to_year_fraction_simple)(value)

    @staticmethod
    def eval(value):
        return mjd2000_to_year_fraction_simple(value)

    @staticmethod
    def _assert(tested, expected):
        assert_allclose(tested, expected, rtol=1e-14, atol=1e-11)

    def test_mjd2000_to_year_fraction_far_range(self):
        values = uniform(-730487., 730485., (100, 100))
        self._assert(self.eval(values), self.reference(values))

    def test_mjd2000_to_year_fraction_near_range(self):
        values = uniform(-36524., 36525., (100, 100))
        self._assert(self.eval(values), self.reference(values))

    def test_mjd2000_to_year_fraction_sanity_check(self):
        self._assert(self.eval(0.), 0.0)
        self._assert(self.eval([-365.25, 365.25]), [0.0, 0.0])
        self._assert(self.eval([6757.125, 7487.625]), [0.5, 0.5])


class TestMjd2000ToDecimalYearSimple(TestCase):

    @staticmethod
    def reference(value):
        return vectorize(_mjd2000_to_decimal_year_simple)(value)

    @staticmethod
    def eval(value):
        return mjd2000_to_decimal_year_simple(value)

    @staticmethod
    def _assert(tested, expected):
        assert_allclose(tested, expected, rtol=1e-14, atol=1e-11)

    def test_mjd2000_to_decimal_year_far_range(self):
        values = uniform(-730487., 730485., (100, 100))
        self._assert(self.eval(values), self.reference(values))

    def test_mjd2000_to_decimal_year_near_range(self):
        values = uniform(-36524., 36525., (100, 100))
        self._assert(self.eval(values), self.reference(values))

    def test_mjd2000_to_decimal_year_sanity_check(self):
        self._assert(self.eval(0.), 2000.0)
        self._assert(self.eval([-365.25, 365.25]), [1999.0, 2001])
        self._assert(self.eval([6757.125, 7487.625]), [2018.5, 2020.5])

    def test_mjd2000_to_decimal_year_special_values(self):
        self._assert(self.eval([-inf, inf, nan]), [-inf, inf, nan])


class TestDecimalYearToMjd2000Simple(TestCase):

    @staticmethod
    def reference(value):
        return vectorize(_decimal_year_to_mjd2000_simple)(value)

    @staticmethod
    def eval(value):
        return decimal_year_to_mjd2000_simple(value)

    @staticmethod
    def _assert(tested, expected):
        assert_allclose(tested, expected, rtol=1e-14, atol=6e-8)

    def test_decimal_year_to_mjd2000_far_range(self):
        values = uniform(0.0, 4000.0, (100, 100))
        self._assert(self.eval(values), self.reference(values))

    def test_decimal_year_to_mjd2000_near_range(self):
        values = uniform(1900, 2100.0, (100, 100))
        self._assert(self.eval(values), self.reference(values))

    def test_decimal_year_to_mjd2000_sanity_check(self):
        self._assert(self.eval(2000.), 0.0)
        self._assert(self.eval([1999., 2001.0]), [-365.25, 365.25])
        self._assert(self.eval([2018.5, 2020.5]), [6757.125, 7487.625])

    def test_decimal_year_to_mjd2000_special_values(self):
        self._assert(self.eval([-inf, inf, nan]), [-inf, inf, nan])

#-------------------------------------------------------------------------------
# reference implementation

def _mjd2000_to_year_fraction_simple(mjd2k):
    """ Convert Modified Julian Date 2000 to (Julian) year fraction.
    """
    return (mjd2k / 365.25) % 1


def _mjd2000_to_decimal_year_simple(mjd2k):
    """ Convert Modified Julian Date 2000 to decimal year.
    """
    return 2000.0 + mjd2k / 365.25


def _decimal_year_to_mjd2000_simple(decimal_year):
    """ Covert decimal year to Modified Julian Date 2000.
    """
    return (decimal_year - 2000.0) * 365.25

#-------------------------------------------------------------------------------

if __name__ == "__main__":
    main()
