#-------------------------------------------------------------------------------
#
#  EMM format parser - test
#
# Author: Martin Paces <martin.paces@eox.at>
#
#-------------------------------------------------------------------------------
# Copyright (C) 2018 EOX IT Services GmbH
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies of this Software or works derived from this Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#-------------------------------------------------------------------------------
# pylint: disable=missing-docstring

from unittest import TestCase, main
from numpy.testing import assert_allclose
from numpy import abs as aabs
from eoxmagmod import decimal_year_to_mjd2000
from eoxmagmod.magnetic_model.parser_emm import (
    combine_emm_coefficients, parse_emm_file, EMM_VALIDITY_PERIOD,
)
from eoxmagmod.data import EMM_2010_STATIC, EMM_2010_SECVAR


class TestEMMParser(TestCase):

    @staticmethod
    def parse(filename):
        with open(filename) as file_in:
            return parse_emm_file(file_in)

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

    def _assert_valid_variable(self, data, expected_data):
        self._assert_valid(data, expected_data)
        assert_allclose(
            data["t"], [data["epoch"], data["epoch"] + EMM_VALIDITY_PERIOD]
        )

    def _assert_valid_constant(self, data, expected_data):
        self._assert_valid(data, expected_data)
        assert_allclose(data["t"], [data["epoch"]])

    def test_parse_emm_file_emm2010_static(self):
        data = self.parse(EMM_2010_STATIC)
        self._assert_valid_constant(data, {
            "epoch": 2010.0,
            "degree_min": 1,
            "degree_max": 740,
        })

    def test_parse_emm_file_emm2010_secvar(self):
        data = self.parse(EMM_2010_SECVAR)
        self._assert_valid_constant(data, {
            "epoch": 2010.0,
            "degree_min": 1,
            "degree_max": 16,
        })

    def test_combine_emm_coefficients_emm2010(self):
        data_variable, data_constant = combine_emm_coefficients(
            self.parse(EMM_2010_STATIC), self.parse(EMM_2010_SECVAR)
        )
        self._assert_valid_variable(data_variable, {
            "epoch": 2010.0,
            "degree_min": 1,
            "degree_max": 16,
        })
        self._assert_valid_constant(data_constant, {
            "epoch": 2010.0,
            "degree_min": 17,
            "degree_max": 739,
        })


if __name__ == "__main__":
    main()
