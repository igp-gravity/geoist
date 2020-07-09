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
# pylint: disable=missing-docstring, invalid-name, too-few-public-methods

from unittest import TestCase, main
from itertools import product
from random import random
from math import pi #, sin, cos, sqrt, floor
from numpy import array, empty, sin, cos, sqrt, arctan2, hypot, fabs, copysign
from numpy.testing import assert_allclose
from geoist.magmod._pymm import (
    GEODETIC_ABOVE_WGS84, GEOCENTRIC_SPHERICAL, GEOCENTRIC_CARTESIAN, convert,
)

WGS84_A = 6378.137
WGS84_B = 6356.7523142
WGS84_EPS2 = (1.0 - (WGS84_B/WGS84_A)**2)
DEG2RAD = pi / 180.0
RAD2DEG = 180.0 / pi


class ConvertTestMixIn(object):
    source_coordinate_system = None
    target_coordinate_system = None

    @property
    def coordinates(self):
        raise NotImplementedError

    @classmethod
    def reference_convert(cls, coords):
        raise NotImplementedError

    @classmethod
    def eval_convert(cls, coords):
        return convert(
            coords, cls.source_coordinate_system, cls.target_coordinate_system
        )

    def test_convert(self):
        coords = self.coordinates
        assert_allclose(
            self.eval_convert(coords), self.reference_convert(coords)
        )

    @staticmethod
    def spherical2cartesian(coords):
        latr = coords[..., 0] * DEG2RAD
        lonr = coords[..., 1] * DEG2RAD
        rad = coords[..., 2]
        rad_cos_lat = rad * cos(latr)
        result = empty(coords.shape)
        result[..., 0] = rad_cos_lat * cos(lonr)
        result[..., 1] = rad_cos_lat * sin(lonr)
        result[..., 2] = rad * sin(latr)
        return result

    @staticmethod
    def cartesian2spherical(coords):
        x = coords[..., 0]
        y = coords[..., 1]
        z = coords[..., 2]
        rad_xy = hypot(x, y)
        result = empty(coords.shape)
        result[..., 0] = RAD2DEG * arctan2(z, rad_xy)
        result[..., 1] = RAD2DEG * arctan2(y, x)
        result[..., 2] = hypot(rad_xy, z)
        return result

    @classmethod
    def spherical2geodetic(cls, coords):
        latr = coords[..., 0] * DEG2RAD
        lond = coords[..., 1]
        rad = coords[..., 2]

        result = empty(coords.shape)
        result[..., 1] = lond
        result[..., 0], result[..., 2] = cls.to_geodetic(
            rad*sin(latr), rad*cos(latr)
        )
        return result

    @classmethod
    def cartesian2geodetic(cls, coords):
        x = coords[..., 0]
        y = coords[..., 1]
        z = coords[..., 2]
        result = empty(coords.shape)
        result[..., 1] = RAD2DEG*arctan2(y, x)
        result[..., 0], result[..., 2] = cls.to_geodetic(z, hypot(x, y))
        return result

    @staticmethod
    def to_geodetic(z_coord, hypot_xy):
        """ Get geodetic coordinates calculated by the Ferrarri's solution. """
        ee4 = WGS84_EPS2**2
        pa2 = (hypot_xy / WGS84_A)**2
        zt = (1.0 - WGS84_EPS2) * (z_coord / WGS84_A)**2
        rh = (pa2 + zt - ee4)/6.0
        ss = (0.25*ee4) * zt * pa2
        rh3 = rh**3
        tmp = rh3 + ss + sqrt(ss*(ss+2.0*rh3))
        tt = copysign(fabs(tmp)**(1.0/3.0), tmp)
        uu = rh + tt + rh**2 / tt
        vv = sqrt(uu**2 + ee4*zt)
        ww = (0.5*WGS84_EPS2) * (uu + vv - zt)/vv
        kp = 1.0 + WGS84_EPS2*(sqrt(uu + vv + ww**2) + ww)/(uu + vv)
        zkp = kp * z_coord
        return (
            RAD2DEG*arctan2(zkp, hypot_xy),
            hypot(hypot_xy, zkp)*(1.0/kp - 1.0 + WGS84_EPS2)/WGS84_EPS2
        )

    @classmethod
    def geodetic2spherical(cls, coords):
        latr = coords[..., 0] * DEG2RAD
        lond = coords[..., 1]
        elev = coords[..., 2]

        z_coord, hypot_xy = cls.from_geodetic(latr, elev)

        result = empty(coords.shape)
        result[..., 0] = RAD2DEG * arctan2(z_coord, hypot_xy)
        result[..., 1] = lond
        result[..., 2] = hypot(z_coord, hypot_xy)
        return result

    @classmethod
    def geodetic2cartesian(cls, coords):
        latr = coords[..., 0] * DEG2RAD
        lonr = coords[..., 1] * DEG2RAD
        elev = coords[..., 2]

        z_coord, hypot_xy = cls.from_geodetic(latr, elev)

        result = empty(coords.shape)
        result[..., 0] = hypot_xy * cos(lonr)
        result[..., 1] = hypot_xy * sin(lonr)
        result[..., 2] = z_coord
        return result

    @staticmethod
    def from_geodetic(latitude, elevation):
        sin_lat = sin(latitude)
        cos_lat = cos(latitude)
        rad_curv = WGS84_A / sqrt(1.0 - WGS84_EPS2 * sin_lat**2)
        hypot_xy = (rad_curv + elevation)*cos_lat
        z_coord = (rad_curv*(1.0 - WGS84_EPS2) + elevation)*sin_lat
        return z_coord, hypot_xy

#-------------------------------------------------------------------------------
# sources

class SourceSpherical(object):
    source_coordinate_system = GEOCENTRIC_SPHERICAL

    @property
    def coordinates(self):
        return array([
            (lat, lon, 6371.2*(1.0 + random())) for lat, lon
            in product(range(-90, 91, 5), range(-180, 181, 10))
        ])


class SourceGeodetic(object):
    source_coordinate_system = GEODETIC_ABOVE_WGS84

    @property
    def coordinates(self):
        return array([
            (lat, lon, -50. + 200.0 * random()) for lat, lon
            in product(range(-90, 91, 5), range(-180, 181, 10))
        ])


class SourceCartesian(object):
    source_coordinate_system = GEOCENTRIC_CARTESIAN

    @property
    def coordinates(self):
        return self.spherical2cartesian(array([
            (lat, lon, 6371.2*(1.0 + random())) for lat, lon
            in product(range(-90, 91, 5), range(-180, 181, 10))
        ]))

#-------------------------------------------------------------------------------
# generic targets

class TargetIndetical(object):

    @classmethod
    def reference_convert(cls, coords):
        return coords

#-------------------------------------------------------------------------------
# tests

class TestGeodeticWGS84ToCartesian(TestCase, SourceGeodetic, ConvertTestMixIn):
    source_coordinate_system = GEODETIC_ABOVE_WGS84
    target_coordinate_system = GEOCENTRIC_CARTESIAN

    @classmethod
    def reference_convert(cls, coords):
        return cls.geodetic2cartesian(coords)


class TestGeodeticWGS84ToSpherical(TestCase, SourceGeodetic, ConvertTestMixIn):
    source_coordinate_system = GEODETIC_ABOVE_WGS84
    target_coordinate_system = GEOCENTRIC_SPHERICAL

    @classmethod
    def reference_convert(cls, coords):
        return cls.geodetic2spherical(coords)


class TestGeodeticWGS84ToGeodeticWGS84(TestCase, SourceGeodetic, TargetIndetical, ConvertTestMixIn):
    source_coordinate_system = GEODETIC_ABOVE_WGS84
    target_coordinate_system = GEODETIC_ABOVE_WGS84


class TestSphericalToCartesian(TestCase, SourceSpherical, ConvertTestMixIn):
    target_coordinate_system = GEOCENTRIC_CARTESIAN

    @classmethod
    def reference_convert(cls, coords):
        return cls.spherical2cartesian(coords)


class TestSphericalToSpherical(TestCase, SourceSpherical, TargetIndetical, ConvertTestMixIn):
    target_coordinate_system = GEOCENTRIC_SPHERICAL


class TestSphericalToGeodeticWGS84(TestCase, SourceSpherical, ConvertTestMixIn):
    target_coordinate_system = GEODETIC_ABOVE_WGS84

    @classmethod
    def reference_convert(cls, coords):
        return cls.spherical2geodetic(coords)


class TestCartesianToCartesian(TestCase, SourceCartesian, TargetIndetical, ConvertTestMixIn):
    target_coordinate_system = GEOCENTRIC_CARTESIAN


class TestCartesianToSpherical(TestCase, SourceCartesian, ConvertTestMixIn):
    target_coordinate_system = GEOCENTRIC_SPHERICAL

    @classmethod
    def reference_convert(cls, coords):
        return cls.cartesian2spherical(coords)


class TestCartesianToGeodeticWGS84(TestCase, SourceCartesian, ConvertTestMixIn):
    target_coordinate_system = GEODETIC_ABOVE_WGS84

    @classmethod
    def reference_convert(cls, coords):
        return cls.cartesian2geodetic(coords)


if __name__ == "__main__":
    main()
