#-------------------------------------------------------------------------------
#
#  Coefficient loaders - tests
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
from numpy import inf
from numpy.testing import assert_allclose
from magmod.time_util import (
    decimal_year_to_mjd2000, decimal_year_to_mjd2000_simple,
)
from magmod.data import (
    CHAOS6_CORE_LATEST, CHAOS6_STATIC,
    IGRF11, IGRF12, SIFM, WMM_2015,
    EMM_2010_STATIC, EMM_2010_SECVAR,
    LCS1, MF7,
)
from magmod.magnetic_model.tests.data import (
    SWARM_MMA_SHA_2C_TEST_DATA,
    SWARM_MMA_SHA_2F_TEST_DATA,
    SWARM_MIO_SHA_2_TEST_DATA,
    CHAOS_MMA_TEST_DATA,
)
from magmod.magnetic_model.coefficients import (
    SparseSHCoefficientsTimeDependent,
    SparseSHCoefficientsTimeDependentDecimalYear,
    SparseSHCoefficientsConstant,
    CombinedSHCoefficients,
)
from magmod.magnetic_model.coefficients_mio import (
    SparseSHCoefficientsMIO,
)
from magmod.magnetic_model.loader_shc import (
    load_coeff_shc, load_coeff_shc_combined,
)
from magmod.magnetic_model.loader_igrf import load_coeff_igrf
from magmod.magnetic_model.loader_wmm import load_coeff_wmm
from magmod.magnetic_model.loader_emm import load_coeff_emm
from magmod.magnetic_model.loader_mma import (
    load_coeff_swarm_mma_2c_internal, load_coeff_swarm_mma_2c_external,
    load_coeff_swarm_mma_2f_geo_internal, load_coeff_swarm_mma_2f_geo_external,
    load_coeff_swarm_mma_2f_sm_internal, load_coeff_swarm_mma_2f_sm_external,
)
from magmod.magnetic_model.loader_mio import (
    load_coeff_swarm_mio_internal, load_coeff_swarm_mio_external
)


class CoefficietLoaderTestMixIn(object):
    is_internal = True
    class_ = None
    validity = None
    degree = 0
    min_degree = -1

    @property
    def coeff(self):
        if not hasattr(self, "_coeff"):
            self._coeff = self.load()
        return self._coeff

    @classmethod
    def load(cls):
        raise NotImplementedError

    def test_validity(self):
        assert_allclose(self.coeff.validity, self.validity, rtol=1e-8, atol=1e-8)

    def test_class(self):
        self.assertIsInstance(self.coeff, self.class_)

    def test_degree(self):
        self.assertEqual(self.coeff.degree, self.degree)

    def test_min_degree(self):
        self.assertEqual(self.coeff.min_degree, self.min_degree)

    def test_model_type(self):
        self.assertEqual(self.coeff.is_internal, self.is_internal)


class MIOCoefficietLoaderTestMixIn(CoefficietLoaderTestMixIn):
    class_ = SparseSHCoefficientsMIO
    validity = (-inf, +inf)
    ps_extent = (0, 4, -2, 2)
    lat_ngp = 80.08
    lon_ngp = -72.22
    height = 110.0
    wolf_ratio = 0.014850

    @property
    def params(self):
        if not hasattr(self, "_params"):
            self._params = self.load_params()
        return self._params

    @classmethod
    def load(cls):
        return cls._load()[0]

    @classmethod
    def load_params(cls):
        return cls._load()[1]

    def test_ps_extent(self):
        self.assertEqual(self.coeff.ps_extent, self.ps_extent)

    def test_ngp_coords(self):
        assert_allclose(
            (self.params["lat_NGP"], self.params["lon_NGP"]),
            (self.lat_ngp, self.lon_ngp),
            rtol=1e-8, atol=1e-8
        )

    def test_height(self):
        assert_allclose(self.params["height"], self.height)

    def test_wolf_ratio(self):
        assert_allclose(self.params["wolf_ratio"], self.wolf_ratio)


class ShcTestMixIn(CoefficietLoaderTestMixIn):
    path = None
    kwargs = {}

    @classmethod
    def load(cls):
        return load_coeff_shc(cls.path, **cls.kwargs)


class CombinedShcTestMixIn(CoefficietLoaderTestMixIn):
    class_ = CombinedSHCoefficients
    path_core = None
    path_static = None
    kwargs = {}

    @classmethod
    def load(cls):
        return load_coeff_shc_combined(
            cls.path_core, cls.path_static, **cls.kwargs
        )


class WmmTestMixIn(CoefficietLoaderTestMixIn):
    path = None

    @classmethod
    def load(cls):
        return load_coeff_wmm(cls.path)

#-------------------------------------------------------------------------------

class TestCoeffSIFM(TestCase, ShcTestMixIn):
    class_ = SparseSHCoefficientsTimeDependentDecimalYear
    path = SIFM
    degree = 70
    min_degree = 1
    kwargs = {
        "interpolate_in_decimal_years": True,
    }
    validity = decimal_year_to_mjd2000((2013.4976, 2015.4962))


class TestCoeffIGRF12(TestCase, ShcTestMixIn):
    class_ = SparseSHCoefficientsTimeDependentDecimalYear
    path = IGRF12
    degree = 13
    min_degree = 1
    kwargs = {
        "interpolate_in_decimal_years": True,
    }
    validity = decimal_year_to_mjd2000((1900.0, 2020.0))


class TestCoeffIGRF11(TestCase, CoefficietLoaderTestMixIn):
    class_ = SparseSHCoefficientsTimeDependentDecimalYear
    degree = 13
    min_degree = 1
    validity = decimal_year_to_mjd2000((1900.0, 2015.0))

    @staticmethod
    def load():
        return load_coeff_igrf(IGRF11)


class TestCoeffLCS1(TestCase, ShcTestMixIn):
    class_ = SparseSHCoefficientsConstant
    path = LCS1
    degree = 185
    min_degree = 1
    validity = (-inf, inf)


class TestCoeffMF7(TestCase, ShcTestMixIn):
    class_ = SparseSHCoefficientsConstant
    path = MF7
    degree = 133
    min_degree = 16
    validity = (-inf, inf)

#-------------------------------------------------------------------------------

class TestCoeffCHAOS6Core(TestCase, ShcTestMixIn):
    class_ = SparseSHCoefficientsTimeDependent
    path = CHAOS6_CORE_LATEST
    degree = 20
    min_degree = 1
    validity = decimal_year_to_mjd2000((1997.102, 2019.7002))


class TestCoeffCHAOS6Static(TestCase, ShcTestMixIn):
    class_ = SparseSHCoefficientsConstant
    path = CHAOS6_STATIC
    degree = 110
    min_degree = 21
    validity = (-inf, inf)


class TestCoeffCHAOS6Combined(TestCase, CombinedShcTestMixIn):
    path_core = CHAOS6_CORE_LATEST
    path_static = CHAOS6_STATIC
    degree = 110
    min_degree = 1
    validity = decimal_year_to_mjd2000((1997.1020, 2019.7002))

#-------------------------------------------------------------------------------

class TestCoeffWMM2015(TestCase, WmmTestMixIn):
    class_ = SparseSHCoefficientsTimeDependentDecimalYear
    path = WMM_2015
    degree = 12
    min_degree = 1
    validity = decimal_year_to_mjd2000((2015., 2020.))


class TestCoeffEMM2010(TestCase, CoefficietLoaderTestMixIn):
    class_ = CombinedSHCoefficients
    degree = 739
    min_degree = 1
    validity = decimal_year_to_mjd2000((2010.0, 2015.0))

    @staticmethod
    def load():
        return load_coeff_emm(EMM_2010_STATIC, EMM_2010_SECVAR)

#-------------------------------------------------------------------------------

class TestCoeffMMA2CInternal(TestCase, CoefficietLoaderTestMixIn):
    is_internal = True
    class_ = CombinedSHCoefficients
    degree = 3
    min_degree = 1
    validity = (6179.125, 6209.875)

    @staticmethod
    def load():
        return load_coeff_swarm_mma_2c_internal(SWARM_MMA_SHA_2C_TEST_DATA)


class TestCoeffMMA2CExternal(TestCase, CoefficietLoaderTestMixIn):
    is_internal = False
    class_ = CombinedSHCoefficients
    degree = 2
    min_degree = 1
    validity = (6179.125, 6209.875)

    @staticmethod
    def load():
        return load_coeff_swarm_mma_2c_external(SWARM_MMA_SHA_2C_TEST_DATA)


class TestCoeffChaosMMAInternal(TestCase, CoefficietLoaderTestMixIn):
    is_internal = True
    class_ = CombinedSHCoefficients
    degree = 2
    min_degree = 1
    validity = (6179.00000, 6209.979167)

    @staticmethod
    def load():
        return load_coeff_swarm_mma_2c_internal(CHAOS_MMA_TEST_DATA)


class TestCoeffChaosMMAExternal(TestCase, CoefficietLoaderTestMixIn):
    is_internal = False
    class_ = CombinedSHCoefficients
    degree = 2
    min_degree = 1
    validity = (6179.00000, 6209.979167)

    @staticmethod
    def load():
        return load_coeff_swarm_mma_2c_external(CHAOS_MMA_TEST_DATA)


class TestCoeffMMA2FGeoInternal(TestCase, CoefficietLoaderTestMixIn):
    is_internal = True
    class_ = SparseSHCoefficientsTimeDependent
    degree = 1
    min_degree = 1
    validity = (6179.03125, 6209.96875)

    @staticmethod
    def load():
        return load_coeff_swarm_mma_2f_geo_internal(SWARM_MMA_SHA_2F_TEST_DATA)


class TestCoeffMMA2FGeoExternal(TestCase, CoefficietLoaderTestMixIn):
    is_internal = False
    class_ = SparseSHCoefficientsTimeDependent
    degree = 1
    min_degree = 1
    validity = (6179.03125, 6209.96875)

    @staticmethod
    def load():
        return load_coeff_swarm_mma_2f_geo_external(SWARM_MMA_SHA_2F_TEST_DATA)


class TestCoeffMMA2FSMInternal(TestCase, CoefficietLoaderTestMixIn):
    is_internal = True
    class_ = SparseSHCoefficientsTimeDependent
    degree = 1
    min_degree = 1
    validity = (6179.03125, 6209.96875)

    @staticmethod
    def load():
        return load_coeff_swarm_mma_2f_sm_internal(SWARM_MMA_SHA_2F_TEST_DATA)


class TestCoeffMMA2FSMExternal(TestCase, CoefficietLoaderTestMixIn):
    is_internal = False
    class_ = SparseSHCoefficientsTimeDependent
    degree = 1
    min_degree = 1
    validity = (6179.03125, 6209.96875)

    @staticmethod
    def load():
        return load_coeff_swarm_mma_2f_sm_external(SWARM_MMA_SHA_2F_TEST_DATA)

#-------------------------------------------------------------------------------

class TestCoeffMIOSecondary(TestCase, MIOCoefficietLoaderTestMixIn):
    is_internal = True
    degree = 2
    min_degree = 1
    options = {}

    @classmethod
    def _load(cls):
        return load_coeff_swarm_mio_internal(
            SWARM_MIO_SHA_2_TEST_DATA, **cls.options
        )


class TestCoeffMIOPrimary(TestCase, MIOCoefficietLoaderTestMixIn):
    is_internal = True # default model above the ionosphere
    degree = 2
    min_degree = 1
    options = {}

    @classmethod
    def _load(cls):
        return load_coeff_swarm_mio_external(
            SWARM_MIO_SHA_2_TEST_DATA, **cls.options
        )


class TestCoeffMIOPrimaryBelowIoSph(TestCoeffMIOPrimary):
    is_internal = False
    options = {"above_ionosphere": False}


class TestCoeffMIOPrimaryAboveIoSph(TestCoeffMIOPrimary):
    is_internal = True
    options = {"above_ionosphere": True}

#-------------------------------------------------------------------------------

if __name__ == "__main__":
    main()
