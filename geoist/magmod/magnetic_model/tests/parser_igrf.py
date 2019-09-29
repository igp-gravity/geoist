#-------------------------------------------------------------------------------
#
#  IGRF format parser - test
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
from numpy import abs as aabs
from magmod.magnetic_model.parser_igrf import parse_igrf_file
from magmod.data import IGRF11


class TestIGRFParser(TestCase):

    @staticmethod
    def parse(filename):
        with open(filename) as file_in:
            return parse_igrf_file(file_in)

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

    def test_parse_igrf_file_igrf11(self):
        data = self.parse(IGRF11)
        self._assert_valid(data, {
            "degree_min": 1,
            "degree_max": 13,
        })


if __name__ == "__main__":
    main()
