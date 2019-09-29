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

from math import modf, floor
from numpy import array, vectorize, inf, nan
from numpy.random import uniform
from numpy.testing import assert_allclose
from unittest import TestCase, main
from magmod._pytimeconv import (
    decimal_year_to_mjd2000, mjd2000_to_decimal_year,
    mjd2000_to_year_fraction,
)


class TestMjd2000ToYearFraction(TestCase):

    @staticmethod
    def reference(value):
        return vectorize(_mjd2000_to_year_fraction)(value)

    @staticmethod
    def eval(value):
        return mjd2000_to_year_fraction(value)

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
        self._assert(self.eval([-365.0, 366.0]), [0.0, 0.0])
        self._assert(self.eval([6757.5, 7488.0]), [0.5, 0.5])


class TestMjd2000ToDecimalYear(TestCase):

    @staticmethod
    def reference(value):
        return vectorize(_mjd2000_to_decimal_year)(value)

    @staticmethod
    def eval(value):
        return mjd2000_to_decimal_year(value)

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
        self._assert(self.eval([-365.0, 366.0]), [1999.0, 2001])
        self._assert(self.eval([6757.5, 7488.0]), [2018.5, 2020.5])

    def test_mjd2000_to_decimal_year_special_values(self):
        self._assert(self.eval([-inf, inf, nan]), [-inf, inf, nan])


class TestDecimalYearToMjd2000(TestCase):

    @staticmethod
    def reference(value):
        return vectorize(_decimal_year_to_mjd2000)(value)

    @staticmethod
    def eval(value):
        return decimal_year_to_mjd2000(value)

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
        self._assert(self.eval([1999., 2001.0]), [-365., 366.])
        self._assert(self.eval([2018.5, 2020.5]), [6757.5, 7488.0])

    def test_decimal_year_to_mjd2000_special_values(self):
        self._assert(self.eval([-inf, inf, nan]), [-inf, inf, nan])

#-------------------------------------------------------------------------------
# reference implementation

def _mjd2000_to_year_fraction(mjd2k):
    """ Convert Modified Julian Date 2000 to year fraction.
    """
    year = day2k_to_year(int(floor(mjd2k)))
    return (mjd2k - year_to_day2k(year)) / days_per_year(year)

def _mjd2000_to_decimal_year(mjd2k):
    """ Convert Modified Julian Date 2000 to decimal year.
    """
    year = day2k_to_year(int(floor(mjd2k)))
    return year + (mjd2k - year_to_day2k(year)) / days_per_year(year)


def _decimal_year_to_mjd2000(decimal_year):
    """ Covert decimal year to Modified Julian Date 2000.
    """
    fraction, year = modf(decimal_year)
    year = int(year)
    return year_to_day2k(year) + fraction * days_per_year(year)


def days_per_year(year):
    """ Return integer number of days in the given year."""
    return 365 + is_leap_year(year)


def is_leap_year(year):
    """ Return boolean flag indicating whether the given year is a leap year
    or not.
    """
    if year > 1582:
        return year % 4 == 0 and year % 100 != 0 or year % 400 == 0
    else:
        return year % 4 == 0


def year_to_day2k(year):
    """ Get the date to number of days since 2000-01-01 of the year start. """
    y__ = year + 4799
    day2k = 365*y__ + y__//4 - 2483321
    if year > 1582:
        day2k += y__//400 - y__//100 + 38
    return day2k


def day2k_to_year(day2k):
    """ Convert integer day number since 2000-01-01 to date as (year, month, day)
    tuple.
    """
    d__ = int(day2k) + 2451545
    f__ = d__ + 1401
    if d__ > 2299160:
        f__ += (((4*d__ + 274277)//146097)*3)//4 - 38
    e__ = 4*f__ + 3
    h__ = 5*((e__ % 1461)//4) + 2
    return e__//1461 - 4716 + (13 - (h__//153 + 2)%12)//12

#-------------------------------------------------------------------------------

if __name__ == "__main__":
    main()
