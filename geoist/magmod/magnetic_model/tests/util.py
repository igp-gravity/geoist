#-------------------------------------------------------------------------------
#
#  shared utilities test
#
# Author: Steve Shi Chen <chenshi80@gmail.com>
# 
# Original Author: Martin Paces <martin.paces@eox.at>
#-------------------------------------------------------------------------------
# Copyright (C) 2019 Geoist team
#
#-------------------------------------------------------------------------------


from unittest import TestCase, main
from io import open
from hashlib import md5
from geoist.magmod.magnetic_model.util import parse_file
from geoist.magmod.magnetic_model.tests.data import SWARM_MIO_SHA_2_TEST_DATA


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
