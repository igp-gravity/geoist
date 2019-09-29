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

import sys
try:
    # Python 2
    from itertools import izip as zip
except ImportError:
    pass
from itertools import product
from numpy import array, vectorize
from numpy.random import uniform
from magmod import mjd2000_to_decimal_year
from magmod.quasi_dipole_coordinates import (
    eval_mlt, eval_subsol, eval_qdlatlon_with_base_vectors,
)

EARTH_RADIUS = 6371.2 # km
START_TIME = 0.0    # MJD2000 / 2000-01-01T00:00:00Z
END_TIME = 7305.0   # MJD2000 / 2020-01-01T00:00:00Z


def generate_test_data(file_out):
    """ Generate test dataset. """
    tmp = array([
        (lat, lon) for lat, lon
        in product(range(-90, 91, 5), range(-180, 181, 10))
    ])

    lat, lon = tmp[..., 0], tmp[..., 1]
    rad = uniform(EARTH_RADIUS, 2*EARTH_RADIUS, lat.shape)
    time_mjd2000 = uniform(START_TIME, END_TIME, lat.shape)
    time_decimal_year = vectorize(mjd2000_to_decimal_year)(time_mjd2000)

    qdlat, qdlon, f11, f12, f21, f22, f__ = eval_qdlatlon_with_base_vectors(
        lat, lon, rad, time_decimal_year
    )
    mlt = eval_mlt(qdlon, time_mjd2000)
    sol_lat, sol_lon = eval_subsol(time_mjd2000)

    header = [
        "MJD2000", "DecimalYear", "Latitude", "Longitude", "Radius",
        "QDLatitude", "QDLongitude", "F11", "F12", "F21", "F22", "F",
        "MagneticLocalTime", "SubsolarLatitude", "SubsolarLongitude",
    ]
    records = zip(
        time_mjd2000, time_decimal_year, lat, lon, rad,
        qdlat, qdlon, f11, f12, f21, f22, f__,
        mlt, sol_lat, sol_lon
    )

    file_out.write("\t".join(header) + "\r\n")
    for record in records:
        file_out.write("\t".join("%.14e" % value for value in record) + "\r\n")


if __name__ == "__main__":
    generate_test_data(sys.stdout)
