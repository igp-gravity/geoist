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
from numpy import array
from numpy.testing import assert_allclose
from magmod._pymm import relradpow


class TestRadialPowerSeries(TestCase):

    @staticmethod
    def reference_internal(value, degree):
        """ Reference implementation - internal field. """
        return array([value ** (-i - 2) for i in range(degree + 1)])

    @staticmethod
    def reference_external(value, degree):
        """ Reference implementation - external field. """
        return array([value ** (i - 1) for i in range(degree + 1)])

    def test_invalid_degree(self):
        self.assertRaises(ValueError, relradpow, 1.0, -1, 1.0)

    def test_invalid_radius(self):
        self.assertRaises(ValueError, relradpow, -1.0, 0, 1.0)

    def test_invalid_reference_radius(self):
        self.assertRaises(ValueError, relradpow, 1.0, 0, -1.0)
        self.assertRaises(ValueError, relradpow, 1.0, 0, 0.0)

    def test_relradpow_zero_degree_external(self):
        assert_allclose(relradpow(2.0, 0, 1.0, is_internal=False), [0.5])

    def test_relradpow_zero_degree_internal(self):
        assert_allclose(relradpow(2.0, 0, 1.0, is_internal=True), [0.25])

    def test_relradpow_default(self):
        assert_allclose(
            relradpow(2*6371.2, 10),
            [
                0.25, 0.125, 0.0625, 0.03125, 0.015625, 0.0078125, 0.00390625,
                0.001953125, 0.0009765625, 0.00048828125, 0.000244140625,
            ]
        )
        assert_allclose(
            relradpow(1.1*6371.2, 256),
            self.reference_internal(1.1, 256)
        )

    def test_relradpow_internal(self):
        assert_allclose(
            relradpow(2*6371.2, 10, is_internal=True),
            [
                0.25, 0.125, 0.0625, 0.03125, 0.015625, 0.0078125, 0.00390625,
                0.001953125, 0.0009765625, 0.00048828125, 0.000244140625,
            ]
        )
        assert_allclose(
            relradpow(1.1*6371.2, 256, is_internal=True),
            self.reference_internal(1.1, 256)
        )

    def test_relradpow_external_with_custom_radius(self):
        assert_allclose(
            relradpow(4.0, 10, reference_radius=2.0, is_internal=False),
            [0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512]
        )
        assert_allclose(
            relradpow(1.1, 256, reference_radius=1.0, is_internal=False),
            self.reference_external(1.1, 256)
        )

    def test_relradpow_internal_with_custom_radius(self):
        assert_allclose(
            relradpow(4.0, 10, reference_radius=2.0, is_internal=True),
            [
                0.25, 0.125, 0.0625, 0.03125, 0.015625, 0.0078125, 0.00390625,
                0.001953125, 0.0009765625, 0.00048828125, 0.000244140625,
            ]
        )
        assert_allclose(
            relradpow(1.1, 256, reference_radius=1.0, is_internal=True),
            self.reference_internal(1.1, 256)
        )


if __name__ == "__main__":
    main()
