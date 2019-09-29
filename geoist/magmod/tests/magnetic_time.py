#-------------------------------------------------------------------------------
#
#  Magnetic time calculations - tests
#
# Author: Steve Shi Chen <chenshi80@gmail.com>
# 
# Original Author: Martin Paces <martin.paces@eox.at>
#-------------------------------------------------------------------------------
# Copyright (C) 2019 Geoist team
#
#-------------------------------------------------------------------------------
# pylint: disable=missing-docstring, invalid-name

from unittest import TestCase, main
from numpy import asarray, empty, nditer
from numpy.random import uniform
from numpy.testing import assert_allclose
from magmod.time_util import decimal_year_to_mjd2000
from magmod.solar_position import sunpos
from magmod.dipole_coords import convert_to_dipole
from magmod.magnetic_time import mjd2000_to_magnetic_universal_time


class TestMjd200ToMagneticUniversalTime(TestCase):
    shape = (25, 25)

    @property
    def times(self):
        return uniform(
            decimal_year_to_mjd2000(1990),
            decimal_year_to_mjd2000(2030),
            self.shape
        )

    @property
    def ngp_coords(self):
        return uniform(-90, 90, self.shape), uniform(-90, 90, self.shape)

    @staticmethod
    def eval(times, lats_ngp, lons_ngp, *args, **kwargs):
        return mjd2000_to_magnetic_universal_time(
            times, lats_ngp, lons_ngp, *args, **kwargs
        )

    @classmethod
    def reference(cls, times, lats_ngp, lons_ngp):
        times = asarray(times)
        lats_ngp = asarray(lats_ngp)
        lons_ngp = asarray(lons_ngp)
        results = empty(times.shape)
        iterator = nditer(
            [times, lats_ngp, lons_ngp, results],
            op_flags=[
                ['readonly'], ['readonly'], ['readonly'], ['writeonly'],
            ],
        )
        for time, lat_ngp, lon_ngp, result in iterator:
            result[...] = cls.ref_mjd2000_to_magnetic_universal_time(
                time, lat_ngp, lon_ngp
            )
        return results

    @staticmethod
    def sunpos(time):
        declination, _, hour_angle, _, _ = sunpos(time, 0, 0, rad=0)
        return declination, -hour_angle

    @classmethod
    def ref_mjd2000_to_magnetic_universal_time(cls, time, lat_ngp, lon_ngp):
        lat_sol, lon_sol = cls.sunpos(time)
        _, subsol_dip_lon, _ = convert_to_dipole(
            [lat_sol, lon_sol, 1.0], lat_ngp, lon_ngp
        )
        return (180.0 - subsol_dip_lon) / 15.0

    def test_mjd2000_to_magnetic_universal_time(self):
        times = self.times
        lats_ngp, lons_ngp = self.ngp_coords
        assert_allclose(
            self.eval(times, lats_ngp, lons_ngp),
            self.reference(times, lats_ngp, lons_ngp),
        )

    def test_mjd2000_to_magnetic_universal_time_fixed_pole(self):
        times = self.times
        lats_ngp, lons_ngp = 80.08, -72.22
        assert_allclose(
            self.eval(times, lats_ngp, lons_ngp),
            self.reference(times, lats_ngp, lons_ngp),
        )

    def test_mjd2000_to_magnetic_universal_time_with_extra_subsol_coords(self):
        times = self.times
        lats_ngp, lons_ngp = self.ngp_coords
        lat_sol, lon_sol = self.sunpos(times)
        assert_allclose(
            self.eval(times, lats_ngp, lons_ngp, lat_sol=lat_sol, lon_sol=lon_sol),
            self.reference(times, lats_ngp, lons_ngp),
        )


if __name__ == "__main__":
    main()
