#-------------------------------------------------------------------------------
#
#  Swarm MMA_SHA_2C product file format parser
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
from numpy import abs as aabs
from spacepy import pycdf
from geoist.magmod.magnetic_model.tests.data import (
    SWARM_MMA_SHA_2C_TEST_DATA, SWARM_MMA_SHA_2F_TEST_DATA,
    CHAOS_MMA_TEST_DATA,
)
from geoist.magmod.magnetic_model.parser_mma import (
    read_swarm_mma_2c_internal, read_swarm_mma_2c_external,
    read_swarm_mma_2f_geo_internal, read_swarm_mma_2f_geo_external,
    read_swarm_mma_2f_sm_internal, read_swarm_mma_2f_sm_external,
)


class SwarmMMAParserMixIn(object):

    @property
    def data(self):
        with pycdf.CDF(self.filename) as cdf:
            return self.parse(cdf)

    def _assert_valid(self, coeff_field, data, expected_data):
        tested_data = {
            key: data[key] for key in expected_data
        }
        self.assertEqual(tested_data, expected_data)
        self.assertEqual(data["t"].size, data[coeff_field].shape[1])
        self.assertEqual(data["nm"].shape[0], data[coeff_field].shape[0])
        self.assertEqual(data["nm"].shape[1], 2)
        self.assertEqual(data["nm"][..., 0].min(), data["degree_min"])
        self.assertEqual(data["nm"][..., 0].max(), data["degree_max"])
        self.assertTrue(aabs(data["nm"][..., 1]).max() <= data["degree_max"])


class TestSwarmMMA2CInternalParser(TestCase, SwarmMMAParserMixIn):
    filename = SWARM_MMA_SHA_2C_TEST_DATA

    @staticmethod
    def parse(cdf):
        return read_swarm_mma_2c_internal(cdf)

    def test_read_swarm_mma_2c_internal(self):
        data = self.data
        self._assert_valid("gh", data[0], {
            "degree_min": 1,
            "degree_max": 1,
        })
        self._assert_valid("gh", data[1], {
            "degree_min": 1,
            "degree_max": 3,
        })


class TestSwarmMMA2CExternalParser(TestCase, SwarmMMAParserMixIn):
    filename = SWARM_MMA_SHA_2C_TEST_DATA

    @staticmethod
    def parse(cdf):
        return read_swarm_mma_2c_external(cdf)

    def test_read_swarm_mma_2c_external(self):
        data = self.data
        self._assert_valid("qs", data[0], {
            "degree_min": 1,
            "degree_max": 1,
        })
        self._assert_valid("qs", data[1], {
            "degree_min": 1,
            "degree_max": 2,
        })


class TestChaosMMAInternalParser(TestCase, SwarmMMAParserMixIn):
    filename = CHAOS_MMA_TEST_DATA

    @staticmethod
    def parse(cdf):
        return read_swarm_mma_2c_internal(cdf)

    def test_read_chaos_mma_internal(self):
        data = self.data
        self._assert_valid("gh", data[0], {
            "degree_min": 1,
            "degree_max": 1,
        })
        self._assert_valid("gh", data[1], {
            "degree_min": 1,
            "degree_max": 2,
        })


class TestChaosMMAExternalParser(TestCase, SwarmMMAParserMixIn):
    filename = CHAOS_MMA_TEST_DATA

    @staticmethod
    def parse(cdf):
        return read_swarm_mma_2c_external(cdf)

    def test_read_chaos_mma_external(self):
        data = self.data
        self._assert_valid("qs", data[0], {
            "degree_min": 1,
            "degree_max": 1,
        })
        self._assert_valid("qs", data[1], {
            "degree_min": 1,
            "degree_max": 2,
        })


class TestSwarmMMA2FGeoInternalParser(TestCase, SwarmMMAParserMixIn):
    filename = SWARM_MMA_SHA_2F_TEST_DATA

    @staticmethod
    def parse(cdf):
        return read_swarm_mma_2f_geo_internal(cdf)

    def test_read_swarm_mma_2f_geo_internal(self):
        self._assert_valid("gh", self.data, {
            "degree_min": 1,
            "degree_max": 1,
        })


class TestSwarmMMA2FGeoExternalParser(TestCase, SwarmMMAParserMixIn):
    filename = SWARM_MMA_SHA_2F_TEST_DATA

    @staticmethod
    def parse(cdf):
        return read_swarm_mma_2f_geo_external(cdf)

    def test_read_swarm_mma_2f_geo_external(self):
        self._assert_valid("qs", self.data, {
            "degree_min": 1,
            "degree_max": 1,
        })


class TestSwarmMMA2FSMInternalParser(TestCase, SwarmMMAParserMixIn):
    filename = SWARM_MMA_SHA_2F_TEST_DATA

    @staticmethod
    def parse(cdf):
        return read_swarm_mma_2f_sm_internal(cdf)

    def test_read_swarm_mma_2f_sm_internal(self):
        self._assert_valid("gh", self.data, {
            "degree_min": 1,
            "degree_max": 1,
        })


class TestSwarmMMA2FSMExternalParser(TestCase, SwarmMMAParserMixIn):
    filename = SWARM_MMA_SHA_2F_TEST_DATA

    @staticmethod
    def parse(cdf):
        return read_swarm_mma_2f_sm_external(cdf)

    def test_read_swarm_mma_2f_sm_external(self):
        self._assert_valid("qs", self.data, {
            "degree_min": 1,
            "degree_max": 1,
        })


if __name__ == "__main__":
    main()
