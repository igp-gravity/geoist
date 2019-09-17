#-------------------------------------------------------------------------------
#
#  Swarm MIO_SHA_2* product file format parser - test
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
from eoxmagmod.magnetic_model.parser_mio import parse_swarm_mio_file
from eoxmagmod.magnetic_model.tests.data import SWARM_MIO_SHA_2_TEST_DATA


class TestSwarmMIOParser(TestCase):

    @staticmethod
    def parse(filename):
        with open(filename) as file_in:
            return parse_swarm_mio_file(file_in)

    def _assert_valid(self, data, expected_data):
        tested_data = {
            key: data[key] for key in expected_data
        }
        self.assertEqual(tested_data, expected_data)
        self.assertEqual(data["gh"].shape[0], data["nm"].shape[0])
        self.assertEqual(data["gh"].shape[1], data["smax"] - data["smin"] + 1)
        self.assertEqual(data["gh"].shape[2], data["pmax"] - data["pmin"] + 1)
        self.assertEqual(data["gh"].shape[3], 2)
        self.assertEqual(data["qs"].shape[0], data["nm"].shape[0])
        self.assertEqual(data["qs"].shape[1], data["smax"] - data["smin"] + 1)
        self.assertEqual(data["qs"].shape[2], data["pmax"] - data["pmin"] + 1)
        self.assertEqual(data["qs"].shape[3], 2)
        self.assertEqual(data["nm"].shape[1], 2)
        self.assertEqual(data["nm"][..., 0].min(), data["degree_min"])
        self.assertEqual(data["nm"][..., 0].max(), data["degree_max"])
        self.assertTrue(aabs(data["nm"][..., 1]).max() <= data["degree_max"])

    def test_parse_swarm_mio_file(self):
        data = self.parse(SWARM_MIO_SHA_2_TEST_DATA)
        self._assert_valid(data, {
            'nmax': 2,
            'mmax': 2,
            'pmin': 0,
            'pmax': 4,
            'smin': -2,
            'smax': 2,
            'theta_NGP': 9.92,
            'phi_NGP': 287.78,
            'lat_NGP': 90 - 9.92,
            'lon_NGP': 287.78 - 360,
            'height': 110.0,
            'wolf_ratio': 0.01485,
        })
        # assert the right order of the Fourier series coefficients
        assert_allclose(data["qs"][0], [
            [
                (1.06992462e-01, 2.78620360e-02),
                (9.29096993e-03, 1.26507647e-02),
                (-3.18043762e-04, -5.00690227e-03),
                (6.84429672e-04, 1.35158938e-03),
                (-3.64334443e-04, -3.66483116e-04),
            ],
            [
                (-7.71206502e-03, 2.35570798e-02),
                (1.77405435e-01, 1.72059684e-01),
                (-2.03209994e-02, 1.04397411e-02),
                (-1.37699395e-03, -5.20757899e-03),
                (9.08831738e-06, 6.00505636e-05),
            ],
            [
                (-1.60998056e+00, 0.00000000e+00),
                (8.14897199e-02, 1.34199720e-01),
                (2.01042603e-02, 1.13266859e-02),
                (9.36382899e-03, 3.89937413e-03),
                (-1.06627005e-03, 2.87669628e-04),
            ],
            [
                (-7.71206502e-03, -2.35570798e-02),
                (1.76095354e-01, 5.92158343e-02),
                (-1.97672602e-02, 7.27189270e-03),
                (6.76187665e-04, -3.72221945e-03),
                (-3.55347333e-04, 6.23110565e-04),
            ],
            [
                (1.06992462e-01, -2.78620360e-02),
                (8.23915923e-04, 3.15713545e-02),
                (5.07104430e-03, 7.14672406e-03),
                (-3.33023952e-03, -5.29365498e-04),
                (3.31154026e-04, -5.32021458e-04),
            ],
        ], atol=1e-14)


if __name__ == "__main__":
    main()
