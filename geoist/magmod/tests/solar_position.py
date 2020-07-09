#-------------------------------------------------------------------------------
#
#  Solar position - test
#
# Author: Steve Shi Chen <chenshi80@gmail.com>
# 
# Original Author: Martin Paces <martin.paces@eox.at>
#-------------------------------------------------------------------------------
# Copyright (C) 2019 Geoist team
#
#-------------------------------------------------------------------------------
#pylint: disable=missing-docstring

from unittest import TestCase, main
from io import open
from numpy import array
from numpy.testing import assert_allclose
from geoist.magmod.quasi_dipole_coordinates import eval_subsol
from geoist.magmod.solar_position import sunpos, sunpos_original
from geoist.magmod.tests.data import SUN_POSITION_TEST_DATA


def load_test_data(filename):
    """ Load test data from a tab-separated values file. """
    def _load_test_data(file_in):
        header = next(file_in).strip().split("\t")
        records = array([
            [float(v) for v in line.strip().split("\t")] for line in file_in
        ])
        return {
            variable: records[..., idx] for idx, variable in enumerate(header)
        }

    with open(filename, encoding="ascii") as file_in:
        return _load_test_data(file_in)


class TestSunPosition(TestCase):
    test_data = load_test_data(SUN_POSITION_TEST_DATA)

    def test_sunpos(self):
        decl, rasc, lha, azimuth, zenith = sunpos(
            self.test_data["MJD2000"],
            self.test_data["Latitude"],
            self.test_data["Longitude"],
            self.test_data["Radius"],
        )
        assert_allclose(decl, self.test_data["Declination"], atol=1e-8)
        assert_allclose(rasc, self.test_data["RightAscension"], atol=1e-8)
        assert_allclose(lha, self.test_data["HourAngle"], atol=1e-8)
        assert_allclose(azimuth, self.test_data["Azimuth"], atol=1e-8)
        assert_allclose(zenith, self.test_data["Zenith"], atol=1e-8)

    def test_sunpos_original(self):
        # Note: the original code does handle correctly only days
        #       between 2000-03-01 and 2400-02-29.
        mask = (
            (self.test_data["MJD2000"] >= 60.0) &
            (self.test_data["MJD2000"] < 146157.0)
        )
        decl, rasc, lha, azimuth, zenith = sunpos_original(
            self.test_data["MJD2000"][mask],
            self.test_data["Latitude"][mask],
            self.test_data["Longitude"][mask],
            self.test_data["Radius"][mask],
            pres=0.0,
        )
        assert_allclose(decl, self.test_data["Declination"][mask], atol=1e-8)
        assert_allclose(rasc, self.test_data["RightAscension"][mask], atol=1e-8)
        assert_allclose(lha, self.test_data["HourAngle"][mask], atol=1e-8)
        assert_allclose(azimuth, self.test_data["Azimuth"][mask], atol=1e-8)
        assert_allclose(zenith, self.test_data["Zenith"][mask], atol=1e-8)


    def test_sunpos_subsol_comparison(self):
        # The two different Solar models are expected to be aligned.
        decl, _, gha, _, _, = sunpos(self.test_data["MJD2000"], 0, 0, 0)
        sollat, sollon = eval_subsol(self.test_data["MJD2000"])
        assert_allclose(decl, sollat, atol=4e-1)
        assert_allclose(-gha, sollon, atol=2e-1)


if __name__ == "__main__":
    main()
