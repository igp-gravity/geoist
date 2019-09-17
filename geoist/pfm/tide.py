# -*- coding: utf-8 -*-
"""
 Name        : tide.py
 Created on  : 2018/11/03 17:00
 Author      : Steve Chen <chenshi@cea-igp.ac.cn>
 Affiliation : Institute of Geophysics, CEA.
 Version     : 0.1.0
 Copyright   : Copyright (C) 2018-2020 GEOIST Development Team. All Rights Reserved.
 License     : Distributed under the MIT License. See LICENSE.txt for more info.
 Github      : https://igp-gravity.github.io/
 Description : code for UTC time, 1.16 don't multiply.
"""

"""Calculate gravitational tide based on Longman 1959."""
#stdlib imports
from collections import namedtuple
from datetime import datetime, timedelta
from math import acos, asin, atan, cos, radians, sin, sqrt
#third party imports
import matplotlib.pyplot as plt
import numpy as np
# local imports

class TideModel():
    """Class to encapsulate the Longman 1959 tide model."""

    def __init__(self):
        """Initialize the model object."""
        self.name = 'Model'
        self.results = namedtuple('results', ['model_time', 'gravity_moon',
                                              'gravity_sun', 'gravity_total'])
        self.results.model_time = []
        self.results.gravity_moon = []
        self.results.gravity_sun = []
        self.results.gravity_total = []

    def calculate_julian_century(self, timestamp):
        """Calculate the julian century and hour.

        Take a datetime object and calculate the decimal Julian century and
        floating point hour. This is in reference to noon on December 31,
        1899 as stated in the Longman paper.

        Parameters
        ----------
        timestamp: datetime
            Time stamp to convert

        Returns
        -------
        float, float
            Julian century and hour

        """
        origin_date = datetime(1899, 12, 31, 12, 00, 00)  # Noon Dec 31, 1899
        dt = timestamp - origin_date
        days = dt.days + dt.seconds / 3600. / 24.
        decimal_julian_century = days / 36525
        julian_hour = (timestamp.hour + timestamp.minute / 60. +
                       timestamp.second / 3600.)
        return decimal_julian_century, julian_hour

    def solve_longman(self, lat, lon, alt, time):
        """Solve the tide model.

        Given the location and datetime object, computes the current
        gravitational tide and associated quantities. Latitude and longitude
        and in the traditional decimal notation, altitude is in meters, time
        is a datetime object.

        Parameters
        ----------
        lat : float
            latitude (in degrees)
        lon : float
            longitude (in degrees)
        alt : float
            altitude (in meters)
        time : datetime
            time at which to solve the model

        Returns
        -------
        float, float, float
            lunar, solar, and total gravitational tides

        """
        T, t0 = self.calculate_julian_century(time)

        if t0 < 0:
            t0 += 24.
        if t0 >= 24:
            t0 -= 24.

        mu = 6.673e-8  # Newton's gravitational constant
        M = 7.3537e25  # Mass of the moon in grams
        S = 1.993e33  # Mass of the sun in grams
        e = 0.05490  # Eccentricity of the moon's orbit
        m = 0.074804  # Ratio of mean motion of the sun to that of the moon
        # Mean distance between the centers of the earth and the moon
        c = 3.84402e10
        # Mean distance between centers of the earth and sun in cm
        c1 = 1.495e13
        h2 = 0.612  # Love parameter
        k2 = 0.303  # Love parameter
        a = 6.378270e8  # Earth's equitorial radius in cm
        i = 0.08979719  # (i) Inclination of the moon's orbit to the ecliptic
        # Inclination of the Earth's equator to the ecliptic 23.452 degrees
        omega = radians(23.452)
        # For some reason his lat/lon is defined with W as + and E as -
        L = -1 * lon
        lamb = radians(lat)  # (lambda) Latitude of point P
        H = alt * 100.  # (H) Altitude above sea-level of point P in cm

        # Lunar Calculations
        # (s) Mean longitude of moon in its orbit reckoned
        # from the referred equinox
        s = (4.72000889397 + 8399.70927456 * T + 3.45575191895e-05 * T * T +
             3.49065850399e-08 * T * T * T)
        # (p) Mean longitude of lunar perigee
        p = (5.83515162814 + 71.0180412089 * T + 0.000180108282532 * T * T +
             1.74532925199e-07 * T * T * T)
        # (h) Mean longitude of the sun
        h = 4.88162798259 + 628.331950894 * T + 5.23598775598e-06 * T * T
        # (N) Longitude of the moon's ascending node in its orbit
        # reckoned from the referred equinox
        N = (4.52360161181 - 33.757146295 * T + 3.6264063347e-05 * T * T +
             3.39369576777e-08 * T * T * T)
        # (I) Inclination of the moon's orbit to the equator
        I = acos(cos(omega)*cos(i) - sin(omega)*sin(i)*cos(N))
        # (nu) Longitude in the celestial equator of its intersection
        # A with the moon's orbit
        nu = asin(sin(i)*sin(N)/sin(I))
        # (t) Hour angle of mean sun measured west-ward from
        # the place of observations
        t = radians(15. * (t0 - 12) - L)

        # (chi) right ascension of meridian of place of observations
        # reckoned from A
        chi = t + h - nu
        # cos(alpha) where alpha is defined in eq. 15 and 16
        cos_alpha = cos(N) * cos(nu) + sin(N) * sin(nu) * cos(omega)
        # sin(alpha) where alpha is defined in eq. 15 and 16
        sin_alpha = sin(omega) * sin(N) / sin(I)
        # (alpha) alpha is defined in eq. 15 and 16
        alpha = 2 * atan(sin_alpha / (1 + cos_alpha))
        # (xi) Longitude in the moon's orbit of its ascending
        # intersection with the celestial equator
        xi = N - alpha

        # (sigma) Mean longitude of moon in radians in its orbit
        # reckoned from A
        sigma = s - xi
        # (l) Longitude of moon in its orbit reckoned from its ascending
        # intersection with the equator
        l = (sigma + 2 * e * sin(s - p) + (5. / 4) * e * e * sin(2 * (s - p)) +
             (15. / 4) * m * e * sin(s - 2 * h + p) + (11. / 8) *
             m * m * sin(2 * (s - h)))

        # Sun
        # (p1) Mean longitude of solar perigee
        p1 = (4.90822941839 + 0.0300025492114 * T + 7.85398163397e-06 *
              T * T + 5.3329504922e-08 * T * T * T)
        # (e1) Eccentricity of the Earth's orbit
        e1 = 0.01675104 - 0.00004180 * T - 0.000000126 * T * T
        # (chi1) right ascension of meridian of place of observations
        # reckoned from the vernal equinox
        chi1 = t + h
        # (l1) Longitude of sun in the ecliptic reckoned from the
        # vernal equinox
        l1 = h + 2 * e1 * sin(h - p1)
        # cosine(theta) Theta represents the zenith angle of the moon
        cos_theta = (sin(lamb) * sin(I) * sin(l) + cos(lamb) * (cos(0.5 * I)**2
                     * cos(l - chi) + sin(0.5 * I)**2 * cos(l + chi)))
        # cosine(phi) Phi represents the zenith angle of the run
        cos_phi = (sin(lamb) * sin(omega) * sin(l1) + cos(lamb) *
                   (cos(0.5 * omega)**2 * cos(l1 - chi1) +
                   sin(0.5 * omega)**2 * cos(l1 + chi1)))

        # Distance
        # (C) Distance parameter, equation 34
        C = sqrt(1. / (1 + 0.006738 * sin(lamb)**2))
        # (r) Distance from point P to the center of the Earth
        r = C * a + H
        # (a') Distance parameter, equation 31
        aprime = 1. / (c * (1 - e * e))
        # (a1') Distance parameter, equation 31
        aprime1 = 1. / (c1 * (1 - e1 * e1))
        # (d) Distance between centers of the Earth and the moon
        d = (1. / ((1. / c) + aprime * e * cos(s - p) + aprime * e * e *
             cos(2 * (s - p)) + (15. / 8) * aprime * m * e * cos(s - 2 * h + p)
             + aprime * m * m * cos(2 * (s - h))))
        # (D) Distance between centers of the Earth and the sun
        D = 1. / ((1. / c1) + aprime1 * e1 * cos(h - p1))

        # (gm) Vertical componet of tidal acceleration due to the moon
        gm = ((mu * M * r / (d * d * d)) * (3 * cos_theta**2 - 1) + (3. / 2) *
              (mu * M * r * r / (d * d * d * d)) *
              (5 * cos_theta**3 - 3 * cos_theta))
        # (gs) Vertical componet of tidal acceleration due to the sun
        gs = mu * S * r / (D * D * D) * (3 * cos_phi**2 - 1)

        love = (1 + h2 - 1.5 * k2)
        g0 = (gm + gs) * 1e3 * love
        return gm * 1e3 * love, gs * 1e3 * love, g0

    def run_model(self):
        """Run the model for a range of times.

        Runs the tidal model beginning at start_time with time steps of
        increment seconds for days.

        """
        self.n_steps = int(24 * self.duration * 3600 / self.increment)

        for i in np.arange(self.n_steps):
            time_at_step = (self.start_time +
                            i * timedelta(seconds=self.increment))
            gm, gs, g = self.solve_longman(self.latitude, self.longitude,
                                           self.altitude, time_at_step)
            self.results.model_time.append(time_at_step)
            self.results.gravity_moon.append(gm)
            self.results.gravity_sun.append(gs)
            self.results.gravity_total.append(g)

    def plot(self):
        """Plot the model results.

        Make a simple plot of the gravitational tide results from the
        model run.
        """
        fig = plt.figure(figsize=(12, 6))
        ax1 = plt.subplot(111)
        ax1.set_xlabel(r'Date', fontsize=18)
        ax1.set_ylabel(r'Anomaly [mGal]', fontsize=18)
        ax1.tick_params(axis='both', which='major', labelsize=16)
        ax1.plot_date(self.results.model_time, self.results.gravity_total,
                      '-k', linewidth=2)
        plt.show()
        return fig, ax1

    def write(self, fname):
        """Write model results to file.

        Write results out of a file for later analysis or reading into another
        method for analysis/correction of data.

        Parameters
        ----------
        fname: string
            name of file to save

        """
        t_string = datetime.strftime(self.start_time, '%Y-%m-%dT%H:%M:%S')
        f = open(fname, 'w')
        f.write('Station latitude: {self.latitude}\n')
        f.write('Station longitude: {self.longitude}\n')
        f.write('Station altitude [m]: {self.altitude}\n')
        f.write('Time Increment [s]: {self.increment}\n')
        f.write('Start Time: {t_string}\n')
        f.write('Duration [days]: {self.duration}\n')
        f.write('\nTime,Lunar,Solar,Total\n')
        f.write('YYYY-MM-DDTHH:MM:SS\tmGal\tmGal\tmGal\n')

        for i in np.arange(self.n_steps):
            t_string = datetime.strftime(self.results.model_time[i],
                                         '%Y-%m-%dT%H:%M:%S')
            f.write('{}\t{}\t{}\t{}\n'.format(t_string,
                                              self.results.gravity_moon[i],
                                              self.results.gravity_sun[i],
                                              self.results.gravity_total[i]))
        f.close()
        
        
if __name__ == "__main__":
    g1=TideModel()
    gdate = datetime(2019, 11, 10, 10, 00, 00)   
    g1.duration = 10
    g1.increment = 60
    g1.start_time = gdate
    g1.latitude = 45.0
    g1.longitude = 105.0
    g1.altitude = 0.0
    g1.run_model()
    g1.plot()
    g1.write('d:/et20191110.txt')
#    gm, gs, g0 =g1.solve_longman(45,105,0,gdate)
#    gdate1 = datetime(1996, 10, 31, 8, 10, 00)     #0.0035347472772194496
#    gm, gs, g00 =g1.solve_longman(22.+44/60,90.5,0.,gdate1)
#    print(g00)