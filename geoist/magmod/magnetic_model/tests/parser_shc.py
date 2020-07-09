#-------------------------------------------------------------------------------
#
#  SHC format parser - test
#
# Author: Steve Shi Chen <chenshi80@gmail.com>
# 
# Original Author: Martin Paces <martin.paces@eox.at>
#-------------------------------------------------------------------------------
# Copyright (C) 2019 Geoist team
#
#-------------------------------------------------------------------------------
# pylint: disable=missing-docstring,invalid-name

from unittest import TestCase, main
from numpy import abs as aabs
from numpy.testing import assert_equal
from geoist.magmod.magnetic_model.parser_shc import parse_shc_file, parse_shc_header
from geoist.magmod.data import (
    CHAOS6_CORE_LATEST, CHAOS6_STATIC,
    IGRF12, SIFM, LCS1, MF7,
)


class TestSHCParser(TestCase):

    @staticmethod
    def parse(filename):
        with open(filename) as file_in:
            return parse_shc_file(file_in)

    @staticmethod
    def parse_header(filename):
        with open(filename) as file_in:
            return parse_shc_header(file_in)

    def _assert_valid(self, data, expected_data):
        tested_data = {
            key: data[key] for key in expected_data
        }
        self.assertEqual(tested_data, expected_data)
        self.assertEqual(data["t"].size, data["gh"].shape[1])
        self.assertEqual(data["nm"].shape[0], data["gh"].shape[0])
        self.assertEqual(data["nm"].shape[1], 2)
        self.assertEqual(data["nm"][..., 0].min(), data["degree_min"])
        self.assertEqual(data["nm"][..., 0].max(), data["degree_max"])
        self.assertTrue(aabs(data["nm"][..., 1]).max() <= data["degree_max"])

    def _test_header(self, filename):
        reference = self.parse(filename)
        tested = self.parse_header(filename)
        assert_equal(reference["t"], tested["t"])
        for remove_key in ["nm", "gh", "t"]:
            reference.pop(remove_key)
        tested.pop("t")
        self.assertEqual(reference, tested)

    def test_parse_shc_file_sifm(self):
        data = self.parse(SIFM)
        self._assert_valid(data, {
            "degree_min": 1,
            "degree_max": 70,
            "spline_order": 2,
            "ntime": 2,
            "nstep": 1,
        })

    def test_parse_shc_header_sifm(self):
        self._test_header(SIFM)

    def test_parse_shc_file_igrf12(self):
        data = self.parse(IGRF12)
        self._assert_valid(data, {
            "degree_min": 1,
            "degree_max": 13,
            "spline_order": 2,
            "ntime": 25,
            "nstep": 1,
        })

    def test_parse_shc_header_igrf12(self):
        self._test_header(IGRF12)

    def test_parse_shc_file_chaos6core_latest(self):
        data = self.parse(CHAOS6_CORE_LATEST)
        self._assert_valid(data, {
            "degree_min": 1,
            "degree_max": 20,
            "spline_order": 6,
            "ntime": 227,
            "nstep": 5,
        })

    def test_parse_shc_header_chaos6core_latest(self):
        self._test_header(CHAOS6_CORE_LATEST)

    def test_parse_shc_file_chaos6static(self):
        data = self.parse(CHAOS6_STATIC)
        self._assert_valid(data, {
            "degree_min": 21,
            "degree_max": 110,
            "spline_order": 1,
            "ntime": 1,
            "nstep": 1,
        })

    def test_parse_shc_header_chaos6static(self):
        self._test_header(CHAOS6_STATIC)

    def test_parse_shc_file_lcs1(self):
        data = self.parse(LCS1)
        self._assert_valid(data, {
            "degree_min": 1,
            "degree_max": 185,
            "spline_order": 1,
            "ntime": 1,
            "nstep": 1,
        })

    def test_parse_shc_file_mf7(self):
        data = self.parse(MF7)
        self._assert_valid(data, {
            "degree_min": 16,
            "degree_max": 133,
            "spline_order": 1,
            "ntime": 1,
            "nstep": 1,
        })


if __name__ == "__main__":
    main()
