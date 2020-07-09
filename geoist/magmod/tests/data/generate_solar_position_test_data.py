#-------------------------------------------------------------------------------
#
#  Solar position - test
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
from numpy import array
from numpy.random import uniform
from geoist.magmod.solar_position import sunpos

EARTH_RADIUS = 6371.2   # km
START_TIME = -10957.0   # MJD2000 / 1970-01-01T00:00:00Z
END_TIME = 10958.0      # MJD2000 / 2030-01-01T00:00:00Z


def generate_test_data(file_out):
    """ Generate test dataset. """
    tmp = array([
        (lat, lon) for lat, lon
        in product(range(-90, 91, 5), range(-180, 181, 10))
    ])
    lat, lon = tmp[..., 0], tmp[..., 1]
    rad = uniform(0*EARTH_RADIUS, 2*EARTH_RADIUS, lat.shape)
    time_mjd2000 = uniform(START_TIME, END_TIME, lat.shape)
    decl, rasc, lha, azim, znth = sunpos(time_mjd2000, lat, lon, rad)

    header = [
        "MJD2000", "Latitude", "Longitude", "Radius",
        "Declination", "RightAscension", "HourAngle", "Azimuth", "Zenith",
    ]
    records = zip(time_mjd2000, lat, lon, rad, decl, rasc, lha, azim, znth)

    file_out.write("\t".join(header) + "\r\n")
    for record in records:
        file_out.write("\t".join("%.14e" % value for value in record) + "\r\n")


if __name__ == "__main__":
    generate_test_data(sys.stdout)
