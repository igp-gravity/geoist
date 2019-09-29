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
# pylint: disable=missing-docstring,no-name-in-module,too-few-public-methods,line-too-long

from unittest import TestCase, main
from itertools import product
from random import random
from numpy import array, empty, nditer
from numpy.testing import assert_allclose
from magmod._pymm import (
    POTENTIAL, GRADIENT, POTENTIAL_AND_GRADIENT,
    GEODETIC_ABOVE_WGS84, GEOCENTRIC_SPHERICAL, GEOCENTRIC_CARTESIAN,
    convert, vrot_sph2geod, vrot_sph2cart,
    relradpow, lonsincos, legendre,
    spharpot, sphargrd, sheval,
)
from magmod.tests.data import sifm
from magmod.tests.data import mma_external


class SphericalHarmonicsMixIn(object):
    options = {}
    scale_potential = 1.0
    scale_gradient = [1.0, 1.0, 1.0]
    source_coordinate_system = None
    target_coordinate_system = None
    is_internal = True
    degree = None
    coef_g = None
    coef_h = None

    @classmethod
    def get_series(cls, degree, latitude, longitude, radius):
        rad_series = relradpow(radius, degree, is_internal=cls.is_internal)
        sin_series, cos_series = lonsincos(longitude, degree)
        p_series, dp_series = legendre(latitude, degree)
        return rad_series, sin_series, cos_series, p_series, dp_series

    @classmethod
    def _spherical_harmonics(cls, latitude, longitude, radius):
        degree = cls.degree
        coef_g = cls.coef_g
        coef_h = cls.coef_h
        rad_series, sin_series, cos_series, p_series, dp_series = cls.get_series(
            degree, latitude, longitude, radius
        )
        potential = spharpot(
            radius, degree, coef_g, coef_h, p_series, rad_series,
            sin_series, cos_series,
        )
        gradient = sphargrd(
            latitude, degree, coef_g, coef_h, p_series, dp_series, rad_series,
            sin_series, cos_series, is_internal=cls.is_internal
        )
        return potential, gradient[0], gradient[1], gradient[2]

    @classmethod
    def _rotate_gradient(cls, vectors, coords):
        if cls.target_coordinate_system == GEOCENTRIC_SPHERICAL:
            return vectors
        elif cls.target_coordinate_system == GEOCENTRIC_CARTESIAN:
            latd = coords[..., 0]
            lond = coords[..., 1]
            return vrot_sph2cart(vectors, latd, lond)
        elif cls.target_coordinate_system == GEODETIC_ABOVE_WGS84:
            dlatd = convert(
                coords, GEOCENTRIC_SPHERICAL, cls.target_coordinate_system
            )[..., 0] - coords[..., 0]
            return vrot_sph2geod(vectors, dlatd)

    @classmethod
    def reference_sheval(cls, coords):
        coords_spherical = convert(
            coords, cls.source_coordinate_system, GEOCENTRIC_SPHERICAL
        )
        potential = empty(coords_spherical.shape[:-1])
        gradient = empty(coords_spherical.shape)

        iterator = nditer(
            [
                potential,
                gradient[..., 0],
                gradient[..., 1],
                gradient[..., 2],
                coords_spherical[..., 0],
                coords_spherical[..., 1],
                coords_spherical[..., 2],
            ],
            op_flags=[
                ['writeonly'], ['writeonly'], ['writeonly'], ['writeonly'],
                ['readonly'], ['readonly'], ['readonly'],
            ],
        )

        for pot, grd_n, grd_e, grd_r, lat, lon, rad in iterator:
            pot[...], grd_n[...], grd_e[...], grd_r[...] = (
                cls._spherical_harmonics(lat, lon, rad)
            )

        gradient = cls._rotate_gradient(gradient, coords_spherical)
        potential *= cls.scale_potential
        for idx, scale in enumerate(cls.scale_gradient):
            gradient[..., idx] *= scale

        return potential, gradient

    @classmethod
    def eval_sheval(cls, coords, mode):
        return sheval(
            coords, mode=mode, is_internal=cls.is_internal,
            degree=cls.degree, coef_g=cls.coef_g, coef_h=cls.coef_h,
            coord_type_in=cls.source_coordinate_system,
            coord_type_out=cls.target_coordinate_system,
            **cls.options
        )

    def test_sheval_potential_and_gradient(self):
        coords = self.coordinates
        potential_ref, gradient_ref = self.reference_sheval(coords)
        potential, gradient = self.eval_sheval(coords, POTENTIAL_AND_GRADIENT)
        assert_allclose(potential, potential_ref, atol=1e-6)
        assert_allclose(gradient, gradient_ref, atol=1e-6)

    def test_sheval_potential(self):
        coords = self.coordinates
        potential_ref, _ = self.reference_sheval(coords)
        potential = self.eval_sheval(coords, POTENTIAL)
        assert_allclose(potential, potential_ref, atol=1e-6)

    def test_sheval_gradient(self):
        coords = self.coordinates
        _, gradient_ref = self.reference_sheval(coords)
        gradient = self.eval_sheval(coords, GRADIENT)
        assert_allclose(gradient, gradient_ref, atol=1e-6)

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
        return convert(array([
            (lat, lon, 6371.2*(1.0 + random())) for lat, lon
            in product(range(-90, 91, 5), range(-180, 181, 10))
        ]), GEOCENTRIC_SPHERICAL, GEOCENTRIC_CARTESIAN)

#-------------------------------------------------------------------------------

class SHTypeInternal(object):
    is_internal = True
    degree = sifm.DEGREE
    coef_g = sifm.COEF_G
    coef_h = sifm.COEF_H


class SHTypeExternal(object):
    is_internal = False
    degree = mma_external.DEGREE
    coef_g = mma_external.COEF_Q
    coef_h = mma_external.COEF_S

#-------------------------------------------------------------------------------

class TestSHEvalCartesian2CartesianInternal(TestCase, SourceCartesian, SHTypeInternal, SphericalHarmonicsMixIn):
    target_coordinate_system = GEOCENTRIC_CARTESIAN

class TestSHEvalCartesian2CartesianExternal(TestCase, SourceCartesian, SHTypeExternal, SphericalHarmonicsMixIn):
    target_coordinate_system = GEOCENTRIC_CARTESIAN


class TestSHEvalCartesian2SphericalInternal(TestCase, SourceCartesian, SHTypeInternal, SphericalHarmonicsMixIn):
    target_coordinate_system = GEOCENTRIC_SPHERICAL

class TestSHEvalCartesian2SphericalExternal(TestCase, SourceCartesian, SHTypeExternal, SphericalHarmonicsMixIn):
    target_coordinate_system = GEOCENTRIC_SPHERICAL


class TestSHEvalCartesian2WGS84Internal(TestCase, SourceCartesian, SHTypeInternal, SphericalHarmonicsMixIn):
    target_coordinate_system = GEODETIC_ABOVE_WGS84

class TestSHEvalCartesian2WGS84External(TestCase, SourceCartesian, SHTypeExternal, SphericalHarmonicsMixIn):
    target_coordinate_system = GEODETIC_ABOVE_WGS84

#-------------------------------------------------------------------------------

class TestSHEvalSpherical2CartesianInternal(TestCase, SourceSpherical, SHTypeInternal, SphericalHarmonicsMixIn):
    target_coordinate_system = GEOCENTRIC_CARTESIAN

class TestSHEvalSpherical2CartesianExternal(TestCase, SourceSpherical, SHTypeExternal, SphericalHarmonicsMixIn):
    target_coordinate_system = GEOCENTRIC_CARTESIAN


class TestSHEvalSpherical2SphericalInternal(TestCase, SourceSpherical, SHTypeInternal, SphericalHarmonicsMixIn):
    target_coordinate_system = GEOCENTRIC_SPHERICAL

class TestSHEvalSpherical2SphericalExternal(TestCase, SourceSpherical, SHTypeExternal, SphericalHarmonicsMixIn):
    target_coordinate_system = GEOCENTRIC_SPHERICAL


class TestSHEvalSpherical2WGS84Internal(TestCase, SourceSpherical, SHTypeInternal, SphericalHarmonicsMixIn):
    target_coordinate_system = GEODETIC_ABOVE_WGS84

class TestSHEvalSpherical2WGS84External(TestCase, SourceSpherical, SHTypeExternal, SphericalHarmonicsMixIn):
    target_coordinate_system = GEODETIC_ABOVE_WGS84

#-------------------------------------------------------------------------------

class TestSHEvalWGS842CartesianInternal(TestCase, SourceGeodetic, SHTypeInternal, SphericalHarmonicsMixIn):
    target_coordinate_system = GEOCENTRIC_CARTESIAN

class TestSHEvalWGS842CartesianExternal(TestCase, SourceGeodetic, SHTypeExternal, SphericalHarmonicsMixIn):
    target_coordinate_system = GEOCENTRIC_CARTESIAN


class TestSHEvalWGS842SphericalInternal(TestCase, SourceGeodetic, SHTypeInternal, SphericalHarmonicsMixIn):
    target_coordinate_system = GEOCENTRIC_SPHERICAL

class TestSHEvalWGS842SphericalExternal(TestCase, SourceGeodetic, SHTypeExternal, SphericalHarmonicsMixIn):
    target_coordinate_system = GEOCENTRIC_SPHERICAL


class TestSHEvalWGS842WGS84Internal(TestCase, SourceGeodetic, SHTypeInternal, SphericalHarmonicsMixIn):
    target_coordinate_system = GEODETIC_ABOVE_WGS84

class TestSHEvalWGS842WGS84External(TestCase, SourceGeodetic, SHTypeExternal, SphericalHarmonicsMixIn):
    target_coordinate_system = GEODETIC_ABOVE_WGS84

#-------------------------------------------------------------------------------

class TestSHEvalCart2CartScaled(TestSHEvalCartesian2CartesianInternal):
    options = {"scale_potential": 2.0, "scale_gradient": [0.5, 1.0, -1.0]}
    scale_potential = 2.0
    scale_gradient = [0.5, 1.0, -1.0]

class TestSHEvalSph2SphScaled(TestSHEvalSpherical2SphericalInternal):
    options = {"scale_gradient": -1.0}
    scale_gradient = [-1.0, -1.0, -1.0]

#-------------------------------------------------------------------------------

if __name__ == "__main__":
    main()
