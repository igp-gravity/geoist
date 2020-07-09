#-------------------------------------------------------------------------------
#
#  Quasi-Dipole coordinates - test
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
from geoist.magmod.quasi_dipole_coordinates import (
    eval_qdlatlon, eval_mlt, eval_subsol,
    eval_qdlatlon_with_base_vectors,
)
from geoist.magmod.tests.data import QUASI_DIPOLE_TEST_DATA


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


class TestQuasiDipoleCoordinates(TestCase):
    test_data = load_test_data(QUASI_DIPOLE_TEST_DATA)

    def test_eval_qdlatlon(self):
        qdlat, qdlon = eval_qdlatlon(
            self.test_data["Latitude"],
            self.test_data["Longitude"],
            self.test_data["Radius"],
            self.test_data["DecimalYear"],
        )
        assert_allclose(qdlat, self.test_data["QDLatitude"], atol=1e-8)
        assert_allclose(qdlon, self.test_data["QDLongitude"], atol=1e-8)

    def test_eval_qdlatlon_with_base_vectors(self):
        qdlat, qdlon, f11, f12, f21, f22, f__ = eval_qdlatlon_with_base_vectors(
            self.test_data["Latitude"],
            self.test_data["Longitude"],
            self.test_data["Radius"],
            self.test_data["DecimalYear"],
        )
        assert_allclose(qdlat, self.test_data["QDLatitude"], atol=1e-8)
        assert_allclose(qdlon, self.test_data["QDLongitude"], atol=1e-8)
        assert_allclose(f11, self.test_data["F11"], atol=1e-8)
        assert_allclose(f12, self.test_data["F12"], atol=1e-8)
        assert_allclose(f21, self.test_data["F21"], atol=1e-8)
        assert_allclose(f22, self.test_data["F22"], atol=1e-8)
        assert_allclose(f__, self.test_data["F"], atol=1e-8)

    def test_eval_mlt(self):
        mlt = eval_mlt(self.test_data["QDLongitude"], self.test_data["MJD2000"])
        assert_allclose(mlt, self.test_data["MagneticLocalTime"], atol=1e-8)

    def test_eval_subsol(self):
        sollat, sollon = eval_subsol(self.test_data["MJD2000"])
        assert_allclose(sollat, self.test_data["SubsolarLatitude"], atol=1e-8)
        assert_allclose(sollon, self.test_data["SubsolarLongitude"], atol=1e-8)


if __name__ == "__main__":
    main()
