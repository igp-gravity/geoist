#-------------------------------------------------------------------------------
#
#  Magnetic Dipole Coordinates - tests
#
# Author: Steve Shi Chen <chenshi80@gmail.com>
# 
# Original Author: Martin Paces <martin.paces@eox.at>
#-------------------------------------------------------------------------------
# Copyright (C) 2019 Geoist team
#
#-------------------------------------------------------------------------------
# pylint: disable=missing-docstring,invalid-name,no-name-in-module

from unittest import TestCase, main
from itertools import product
from math import pi
from numpy import array, asarray, zeros, linspace, meshgrid, sin, cos, dot
from numpy.random import random
from numpy.testing import assert_allclose
from geoist.magmod._pymm import (
    convert, vrot_sph2cart, vrot_cart2sph,
    GEOCENTRIC_SPHERICAL, GEOCENTRIC_CARTESIAN, GEODETIC_ABOVE_WGS84,
)
from geoist.magmod.dipole_coords import (
    get_dipole_rotation_matrix,
    convert_to_dipole,
    vrot_from_dipole,
)

DEG2RAD = pi/180.0


class TestDipoleRotationMatrix(TestCase):

    @staticmethod
    def reference_rotation_matrix(latitude, longitude):
        sin_lat, cos_lat = sin(DEG2RAD*latitude), cos(DEG2RAD*latitude)
        sin_lon, cos_lon = sin(DEG2RAD*longitude), cos(DEG2RAD*longitude)
        matrix = dot(
            # rotate around azimuth axis by -longitude
            array([
                [cos_lon, -sin_lon, 0],
                [sin_lon, cos_lon, 0],
                [0, 0, 1],
            ]),
            # rotate around elevation axis by 90dg - latitude
            array([
                [sin_lat, 0, cos_lat],
                [0, 1, 0],
                [-cos_lat, 0, sin_lat],
            ])
        )
        return matrix

    @staticmethod
    def eval_rotation_matrix(latitude, longitude):
        return get_dipole_rotation_matrix(latitude, longitude)

    def test_rotation_matrix(self):
        coords = [
            (lat, lon) for lat, lon
            in product(range(-90, 91, 5), range(-180, 181, 10))
        ]

        for lat, lon in coords:
            matrix = self.eval_rotation_matrix(lat, lon)
            assert_allclose(
                matrix,
                self.reference_rotation_matrix(lat, lon),
                atol=1e-14
            )
            assert_allclose(
                dot(matrix.transpose(), matrix),
                [(1, 0, 0), (0, 1, 0), (0, 0, 1)],
                atol=1e-14
            )


class TestConvertToDipoleCoordinates(TestCase):

    @staticmethod
    def reference_convert_to_dipole(coords, latitude, longitude):
        rotation_matrix = get_dipole_rotation_matrix(latitude, longitude)
        coords = convert(coords, GEOCENTRIC_SPHERICAL, GEOCENTRIC_CARTESIAN)
        coords = dot(coords, rotation_matrix)
        return coords

    @staticmethod
    def eval_convert_to_dipole(coords, latitude, longitude):
        # to avoid pole longitude ambiguity compare Cartesian coordinates
        return convert(
            convert_to_dipole(coords, latitude, longitude),
            GEOCENTRIC_SPHERICAL, GEOCENTRIC_CARTESIAN
        )

    @property
    def coordinates(self):
        return array([
            (lat, lon, 6371.2*(1.0 + random())) for lat, lon
            in product(range(-90, 91, 5), range(-180, 181, 10))
        ])


    def test_convert_to_dipole(self):
        north_pole_coords = [
            (lat, lon) for lat, lon
            in product(range(-90, 91, 10), range(-180, 181, 20))
        ]
        for lat, lon in north_pole_coords:
            coords = self.coordinates
            assert_allclose(
                self.eval_convert_to_dipole(coords, lat, lon),
                self.reference_convert_to_dipole(coords, lat, lon),
                atol=1e-8
            )

    def test_convert_to_dipole_sanity_check(self):
        assert_allclose(
            self.eval_convert_to_dipole([
                (80, -170, 1.0),
                (-80, 10, 1.0),
                (-10, -170, 1.0),
                (10, 10, 1.0),
                (0, -80, 1.0),
                (0, 100, 1.0),
            ], 80, -170),
            [
                (0, 0, 1),
                (0, 0, -1),
                (1, 0, 0),
                (-1, 0, 0),
                (0, 1, 0),
                (0, -1, 0),
            ],
            atol=1e-12
        )


class VRotFromDipoleMixIn(object):
    target_coords_type = None
    shape = (37, 37)

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
        coords[..., 2] = 6371.2
        return coords

    @classmethod
    def reference_vrot_from_dipole(cls, vectors, coords, latitude, longitude):
        coords = asarray(coords)
        rotation_matrix = get_dipole_rotation_matrix(latitude, longitude)
        coords_dipole = convert_to_dipole(coords, latitude, longitude)
        lat_dipole = coords_dipole[..., 0]
        lon_dipole = coords_dipole[..., 1]
        vectors = vrot_sph2cart(vectors, lat_dipole, lon_dipole)
        vectors = dot(vectors, rotation_matrix.transpose())
        if cls.target_coords_type != GEOCENTRIC_CARTESIAN:
            coords_out = convert(
                coords, GEOCENTRIC_SPHERICAL, cls.target_coords_type
            )
            lat_out = coords_out[..., 0]
            lon_out = coords_out[..., 1]
            vectors = vrot_cart2sph(vectors, lat_out, lon_out)
        return vectors

    @classmethod
    def eval_vrot_from_dipole(cls, vectors, coords, latitude, longitude):
        coords = asarray(coords)
        coords_dipole = convert_to_dipole(coords, latitude, longitude)
        lat_dipole = coords_dipole[..., 0]
        lon_dipole = coords_dipole[..., 1]

        if cls.target_coords_type != GEOCENTRIC_CARTESIAN:
            coords_out = convert(
                coords, GEOCENTRIC_SPHERICAL, cls.target_coords_type
            )
            lat_out = coords_out[..., 0]
            lon_out = coords_out[..., 1]
        else:
            lat_out, lon_out = None, None

        return vrot_from_dipole(
            vectors, latitude, longitude, lat_dipole, lon_dipole,
            lat_out, lon_out, cls.target_coords_type
        )

    def test_vrot_dipole2spherical(self):
        north_pole_coords = [
            (lat, lon) for lat, lon
            in product(range(-90, 91, 10), range(-180, 181, 20))
        ]
        for lat, lon in north_pole_coords:
            coords = self.coords
            vects = self.vectors
            assert_allclose(
                self.eval_vrot_from_dipole(vects, coords, lat, lon),
                self.reference_vrot_from_dipole(vects, coords, lat, lon),
                atol=1e-12
            )


class TestVRotDipoleToCartesian(TestCase, VRotFromDipoleMixIn):
    target_coords_type = GEOCENTRIC_CARTESIAN

    def test_vrot_dipole_to_cartesian_sanity_check(self):
        lat_ngp, lon_ngp = 80, -170
        vectors = array([
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1),
        ])
        sin10, cos10 = sin(DEG2RAD*10), cos(DEG2RAD*10)
        input_output_pairs = [
            ((-10, -170), (0, 0), [
                (-sin10*cos10, -sin10**2, cos10),
                (sin10, -cos10, 0),
                (-cos10**2, -sin10*cos10, -sin10),
            ]),
            ((10, 10), (0, 180), [
                (-sin10*cos10, -sin10**2, cos10),
                (-sin10, cos10, 0),
                (cos10**2, sin10*cos10, sin10),
            ]),
            ((0, -80), (0, 90), [
                (-sin10*cos10, -sin10**2, cos10),
                (cos10**2, sin10*cos10, sin10),
                (sin10, -cos10, 0),
            ]),
            ((0, 100), (0, -90), [
                (-sin10*cos10, -sin10**2, cos10),
                (-cos10**2, -sin10*cos10, -sin10),
                (-sin10, cos10, 0),
            ]),
            ((80, -170), (90, 0), [
                (cos10**2, sin10*cos10, sin10),
                (sin10, -cos10, 0),
                (-sin10*cos10, -sin10**2, cos10),
            ]),
            ((-80, 10), (-90, 0), [
                (-cos10**2, -sin10*cos10, -sin10),
                (sin10, -cos10, 0),
                (sin10*cos10, sin10**2, -cos10),
            ]),
        ]
        for (lat_sph, lon_sph), (lat_dip, lon_dip), expected in input_output_pairs:
            assert_allclose(
                vrot_from_dipole(
                    vectors, lat_ngp, lon_ngp, lat_dip, lon_dip,
                    lat_sph, lon_sph, self.target_coords_type,
                ), expected, atol=1e-12
            )


class VRotFromDipoleToSphericalMixIn(object):
    target_coords_type = GEOCENTRIC_SPHERICAL

    def test_vrot_dipole_to_spherical_sanity_check(self):
        lat_ngp, lon_ngp = 80, -170
        vectors = array([
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1),
        ])
        sin10, cos10 = sin(DEG2RAD*10), cos(DEG2RAD*10)
        input_output_pairs = [
            ((-10, -170), (0, 0), [(1, 0, 0), (0, 1, 0), (0, 0, 1),]),
            ((10, 10), (0, 180), [(1, 0, 0), (0, 1, 0), (0, 0, 1),]),
            (
                (0, -80), (0, 90),
                [(cos10, -sin10, 0), (sin10, cos10, 0), (0, 0, 1)]
            ),
            (
                (0, 100), (0, -90),
                [(cos10, sin10, 0), (-sin10, cos10, 0), (0, 0, 1)]
            ),
            ((80, -170), (90, 0), [(1, 0, 0), (0, 1, 0), (0, 0, 1),]),
            ((-80, 10), (-90, 0), [(-1, 0, 0), (0, -1, 0), (0, 0, 1),]),
        ]
        for (lat_sph, lon_sph), (lat_dip, lon_dip), expected in input_output_pairs:
            assert_allclose(
                vrot_from_dipole(
                    vectors, lat_ngp, lon_ngp, lat_dip, lon_dip,
                    lat_sph, lon_sph, self.target_coords_type,
                ), expected, atol=1e-12
            )


class TestVRotDipoleToSpherical(TestCase, VRotFromDipoleToSphericalMixIn):
    target_coords_type = GEOCENTRIC_SPHERICAL


class TestVRotDipoleToWGS84(TestCase, VRotFromDipoleToSphericalMixIn):
    target_coords_type = GEODETIC_ABOVE_WGS84


if __name__ == "__main__":
    main()
