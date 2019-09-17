#-------------------------------------------------------------------------------
#
#  shared utilities test
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

from unittest import TestCase, main
from io import open
from hashlib import md5
from eoxmagmod.magnetic_model.util import parse_file
from eoxmagmod.magnetic_model.tests.data import SWARM_MIO_SHA_2_TEST_DATA


class TestUtil(TestCase):

    @staticmethod
    def _check_sum(file_in):
        md5sum = md5()
        md5sum.update(file_in.read().encode("UTF-8"))
        return md5sum.hexdigest()

    def test_parse_file(self):
        filename = SWARM_MIO_SHA_2_TEST_DATA
        data_file_name = parse_file(self._check_sum, filename)
        with open(filename, encoding="ascii") as file_in:
            data_file_object = parse_file(self._check_sum, file_in)
        self.assertEqual(data_file_name, data_file_object)


if __name__ == "__main__":
    main()
