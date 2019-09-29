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
from math import pi
from numpy import array, zeros, empty, meshgrid, linspace, sin, cos
from numpy.random import random
from numpy.testing import assert_allclose
from magmod._pymm import (
    vrot_sph2geod, vrot_sph2cart, vrot_cart2sph,
    convert, GEODETIC_ABOVE_WGS84, GEOCENTRIC_SPHERICAL,
)

DEG2RAD = pi/180.0


class VectorRotationMixIn(object):
    shape = (19, 19)

    @property
    def vectors(self):
        return 2.0*random(self.shape + (3,)) - 1.0

    @property
    def coords(self):
        coords = zeros(self.shape + (3,))
        coords[..., 1], coords[..., 0] = meshgrid(
            linspace(-180, 180, self.shape[1]),
            linspace(-90, 90, self.shape[0])
        )
        return coords

    @staticmethod
    def rotate2d(x_coords, y_coords, angles):
        sin_angles, cos_angles = sin(angles), cos(angles)
        return (
            x_coords*cos_angles - y_coords*sin_angles,
            x_coords*sin_angles + y_coords*cos_angles,
        )

    def test_vrot(self):
        vectors = self.vectors
        coords = self.coords
        assert_allclose(
            self.eval_vrot(vectors, coords),
            self.reference_vrot(vectors, coords)
        )

    def test_vrot_sanity_check(self):
        vectors = array([
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 0.0, 1.0),
        ])
        assert_allclose(vrot_cart2sph(vectors, 90, 90), array([
            (0.0, -1.0, 0.0),
            (-1.0, 0.0, 0.0),
            (0.0, 0.0, 1.0),
        ]), atol=1e-14)
        assert_allclose(vrot_cart2sph(vectors, 90, 0), array([
            (-1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 0.0, 1.0),
        ]), atol=1e-14)
        assert_allclose(vrot_cart2sph(vectors, 0, 90), array([
            (0.0, -1.0, 0.0),
            (0.0, 0.0, 1.0),
            (1.0, 0.0, 0.0),
        ]), atol=1e-14)


class TestVectorRotationSphericalToCartesian(TestCase, VectorRotationMixIn):

    @classmethod
    def reference_vrot(cls, vectors, coords):
        latr = coords[..., 0]*DEG2RAD
        lonr = coords[..., 1]*DEG2RAD
        v_n = vectors[..., 0]
        v_e = vectors[..., 1]
        v_r = vectors[..., 2]
        tmp, v_z = cls.rotate2d(v_r, v_n, latr)
        v_x, v_y = cls.rotate2d(tmp, v_e, lonr)
        result = empty(vectors.shape)
        result[..., 0] = v_x
        result[..., 1] = v_y
        result[..., 2] = v_z
        return result

    @staticmethod
    def eval_vrot(vectors, coords):
        latd = coords[..., 0]
        lond = coords[..., 1]
        return vrot_sph2cart(vectors, latd, lond)

    def test_vrot_sanity_check(self):
        vectors = array([
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 0.0, 1.0),
        ])
        assert_allclose(vrot_sph2cart(vectors, 90, 90), array([
            (0.0, -1.0, 0.0),
            (-1.0, 0.0, 0.0),
            (0.0, 0.0, 1.0),
        ]), atol=1e-14)
        assert_allclose(vrot_sph2cart(vectors, 90, 0), array([
            (-1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 0.0, 1.0),
        ]), atol=1e-14)
        assert_allclose(vrot_sph2cart(vectors, 0, 90), array([
            (0.0, 0.0, 1.0),
            (-1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
        ]), atol=1e-14)


class TestVectorRotationCaresianToSpherical(TestCase, VectorRotationMixIn):

    @classmethod
    def reference_vrot(cls, vectors, coords):
        latr = coords[..., 0]*DEG2RAD
        lonr = coords[..., 1]*DEG2RAD
        v_x = vectors[..., 0]
        v_y = vectors[..., 1]
        v_z = vectors[..., 2]
        tmp, v_e = cls.rotate2d(v_x, v_y, -lonr)
        v_r, v_n = cls.rotate2d(tmp, v_z, -latr)
        result = empty(vectors.shape)
        result[..., 0] = v_n
        result[..., 1] = v_e
        result[..., 2] = v_r
        return result

    @staticmethod
    def eval_vrot(vectors, coords):
        latd = coords[..., 0]
        lond = coords[..., 1]
        return vrot_cart2sph(vectors, latd, lond)


class TestVectorRotationSphericalToGeodetic(TestCase, VectorRotationMixIn):

    @staticmethod
    def difflat_spherical_to_geodetic(latitude):
        coords = zeros(latitude.shape + (3,))
        coords[..., 0] = latitude
        coords[..., 2] = 6371.2
        return convert(
            coords, GEOCENTRIC_SPHERICAL, GEODETIC_ABOVE_WGS84
        )[..., 0] - latitude

    @classmethod
    def reference_vrot(cls, vectors, coords):
        dlat_s2g = cls.difflat_spherical_to_geodetic(coords[..., 0])
        v_sn = vectors[..., 0]
        v_se = vectors[..., 1]
        v_sr = vectors[..., 2]
        v_gr, v_gn = cls.rotate2d(v_sr, v_sn, -DEG2RAD*dlat_s2g)
        result = empty(vectors.shape)
        result[..., 0] = v_gn
        result[..., 1] = v_se
        result[..., 2] = v_gr
        return result

    @classmethod
    def eval_vrot(cls, vectors, coords):
        dlat_s2g = cls.difflat_spherical_to_geodetic(coords[..., 0])
        return vrot_sph2geod(vectors, dlat_s2g)

    def test_vrot_sanity_check(self):
        vectors = array([
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 0.0, 1.0),
        ])
        expected_result = array([
            (0.0, 0.0, -1.0),
            (0.0, 1.0, 0.0),
            (1.0, 0.0, 0.0),
        ])
        assert_allclose(vrot_sph2geod(vectors, -90), expected_result, atol=1e-14)


class TestVectorRotationGeodeticToSpherical(TestCase, VectorRotationMixIn):

    @staticmethod
    def difflat_geodetic_to_spherical(latitude):
        coords = zeros(latitude.shape + (3,))
        coords[..., 0] = latitude
        return convert(
            coords, GEODETIC_ABOVE_WGS84, GEOCENTRIC_SPHERICAL
        )[..., 0] - latitude

    @classmethod
    def reference_vrot(cls, vectors, coords):
        dlat_g2s = cls.difflat_geodetic_to_spherical(coords[..., 0])
        v_gn = vectors[..., 0]
        v_ge = vectors[..., 1]
        v_gr = vectors[..., 2]
        v_sr, v_sn = cls.rotate2d(v_gr, v_gn, -DEG2RAD*dlat_g2s)
        result = empty(vectors.shape)
        result[..., 0] = v_sn
        result[..., 1] = v_ge
        result[..., 2] = v_sr
        return result

    @classmethod
    def eval_vrot(cls, vectors, coords):
        dlat_g2s = cls.difflat_geodetic_to_spherical(coords[..., 0])
        return vrot_sph2geod(vectors, dlat_g2s)

    def test_vrot_sanity_check(self):
        vectors = array([
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 0.0, 1.0),
        ])
        expected_result = array([
            (0.0, 0.0, 1.0),
            (0.0, 1.0, 0.0),
            (-1.0, 0.0, 0.0),
        ])
        assert_allclose(vrot_sph2geod(vectors, 90), expected_result, atol=1e-14)


if __name__ == "__main__":
    main()
