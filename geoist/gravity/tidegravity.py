# -*- coding: utf-8 -*-

from datetime import datetime, timezone
from math import sin, radians
from typing import Union
import numpy as np
import pandas as pd


__all__ = ['calculate_julian_century', 'solve_longman_tide', 'solve_longman_tide_scalar',
           'solve_tide_df', 'solve_point_corr']


"""
Longman Earth Tide Calculator - Adopted from John Leeman's implementation of 
I. M. Longman's earth tide calculations (see references below) 

Parts of this program (c) 2017 John Leeman
Any other modifications from the base LongmanTide (c) 2018 Zachery Brady

Licensed under the MIT License, see LICENSE file for full text

References
----------
I.M. Longman "Forumlas for Computing the Tidal Accelerations Due to the Moon 
and the Sun" Journal of Geophysical Research, vol. 64, no. 12, 1959, 
pp. 2351-2355

P. Schureman "Manual of harmonic analysis and prediction of tides" U.S. Coast and Geodetic Survey, 1958

John Leeman's GitHub page for the original implementation of this library: https://github.com/jrleeman/LongmanTide

Notes
-----
Unicode greek symbols are used to more clearly name variables based on the equations in Longman's paper.
This is simply a style decision and may not reflect Python best practices.

"""

# Constants Definitions #
μ = 6.673e-8  # Newton's gravitational constant in cgs units (Verify this, should be 6.674e-8?)
# μ = 6.67408e-8  # Newton's gravitational constant in cgs units (from NIST.gov)
M = 7.3537e25  # Mass of the moon in grams
S = 1.993e33  # Mass of the sun in grams
e = 0.05490  # Eccentricity of the moon's orbit
m = 0.074804  # Ratio of mean motion of the sun to that of the moon
c = 3.84402e10  # Mean distance between the centers of the earth and the moon in cm
c1 = 1.495e13  # Mean distance between centers of the earth and sun in cm
h2 = 0.612  # Love parameter  # See: https://en.wikipedia.org/wiki/Love_number
k2 = 0.303  # Love parameter  # See Also: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4599447/
love = 1 + h2 - 1.5 * k2
a = 6.378270e8  # Earth's equitorial radius in cm
i = 0.08979719  # (i) Inclination of the moon's orbit to the ecliptic
ω = radians(23.452)  # Inclination of the Earth's equator to the ecliptic
origin_date = datetime(1899, 12, 31, 12, 00, 00)  # Noon Dec 31, 1899
# End Constants Definitions #


def calculate_julian_century(dates: Union[np.ndarray, pd.DatetimeIndex]):
    """
    Calculate the decimal Julian century and floating point hour, as referenced from
    1899, December 31 at 12:00:00
    This function accepts either a numpy ndarray or pandas DatetimeIndex as input,
    and returns a 2-tuple of the corresponding Julian century decimals and floating
    point hours.
    Note: All date/times should be supplied as UTC

    Parameters
    ----------
    dates : Union[np.ndarray, pd.DatetimeIndex]
        1-dimensional array of DateTime objects to convert into
        century/hours arrays

    Returns
    -------
    2-tuple of:
        T : np.ndarray
            Number of Julian centuries (36525 days) from GMT Noon on
            December 31, 1899
        t0 : np.ndarray
            Greenwich civil dates measured in hours

    """
    if isinstance(dates, np.ndarray):
        delta = dates - origin_date
        days = np.array([x.days + x.seconds / 3600. / 24. for x in delta])
        t0 = np.array([x.hour + x.minute / 60. + x.second / 3600. for x in
                       dates])
        return days / 36525, t0
    elif isinstance(dates, pd.DatetimeIndex):
        delta = dates - origin_date
        days = delta.days + delta.seconds / 3600. / 24.
        t0 = dates.hour + dates.minute / 60. + dates.second / 3600.
        return days / 36525, t0.values


def solve_longman_tide(lat: np.ndarray, lon: np.ndarray, alt: np.ndarray, time: np.ndarray):
    """
    Find the total gravity correction due to the Sun/Moon for the given
    latitude, longitude, altitude, and time - supplied as numpy 1-d arrays.
    Corrections are calculated for each corresponding set of data in the
    supplied arrays, all arrays must be of the same shape (length and dimension)

    This function returns the Lunar, Solar, and Total gravity corrections as
    a 3-tuple of numpy 1-d arrays.

    Parameters
    ----------
    lat : np.ndarray
        1-dimensional array of float values denoting latitudes
    lon : np.ndarray
        1-dimensional array of float values denoting longitudes
    alt : np.ndarray
        1-dimensional array of float values denoting altitude in meters
    time : np.ndarray
        1-dimensional array of DateTime objects denoting the time series

    Returns
    -------
    3-Tuple
        gMoon (gm): Vertical component of tidal acceleration due to the moon
        gSun (gs): Vertical component of tidal acceleration due to the sun
        gTotal (g0): Total vertical component of tidal acceleration (moon + sun)

    """
    assert lat.shape == lon.shape == alt.shape == time.shape

    T, t0 = calculate_julian_century(time)
    T2 = T ** 2
    T3 = T ** 3

    # t0 must be an ndarray, not a pandas DatetimeIndex
    t0[t0 < 0] += 24
    t0[t0 >= 24] -= 24

    # longitude is defined with West positive (+) in the Longman paper
    φ = -lon
    λ = np.radians(lat)  # λ Latitude of point P
    cosλ = np.cos(λ)
    sinλ = np.sin(λ)
    H = alt * 100  # height above sea-level of point P in centimeters (cm)

    # Lunar Calculations

    # s Mean longitude of moon in its orbit reckoned from the referred equinox
    # Constants from Bartels [1957 pp. 747] eq (10')
    # 270°26'11.72" + (1336 rev. + 1,108,406.05")T + 7.128" * T2 + 0.0072" * T3
    # Where rev. is revolutions expressed as radians: 1336 rev. = 1336 * 2 * pi = 8394.335570392
    s = 4.72000889397 + 8399.70927456 * T + 3.45575191895e-05 * T2 + 3.49065850399e-08 * T3
    # p Mean longitude of lunar perigee
    # constants from Bartels [1957] eq (11')
    # p = 334° 19' 46.42" + (11 rev. + 392,522.51") T - 37.15" * T2 - 0.036" T3
    p = 5.83515162814 + 71.0180412089 * T + 0.000180108282532 * T2 + 1.74532925199e-07 * T3
    # (h) Mean longitude of the sun
    h = 4.88162798259 + 628.331950894 * T + 5.23598775598e-06 * T2
    # (N) Longitude of the moon's ascending node in its orbit reckoned from the referred equinox
    N = 4.52360161181 - 33.757146295 * T + 3.6264063347e-05 * T2 + 3.39369576777e-08 * T3
    cosN = np.cos(N)
    sinN = np.sin(N)

    # I (uppercase i) Inclination of the moon's orbit to the equator
    I = np.arccos(np.cos(ω) * np.cos(i) - np.sin(ω) * np.sin(i) * cosN)
    # ν (nu) Longitude in the celestial equator of its intersection A with the moon's orbit
    ν = np.arcsin(np.sin(i) * sinN / np.sin(I))
    # t Hour angle of mean sun measured west-ward from the place of observations
    t = np.radians(15. * (t0 - 12) - φ)

    # χ (chi) right ascension of meridian of place of observations reckoned from A
    χ = t + h - ν
    # cos α (alpha) where α is defined in eq. 15 and 16
    cos_α = cosN * np.cos(ν) + sinN * np.sin(ν) * np.cos(ω)
    # sin α (alpha) where α is defined in eq. 15 and 16
    sin_α = sin(ω) * sinN / np.sin(I)
    # (α) α is defined in eq. 15 and 16
    α = 2 * np.arctan(sin_α / (1 + cos_α))
    # ξ (xi) Longitude in the moon's orbit of its ascending intersection with the celestial equator
    ξ = N - α

    # σ (sigma) Mean longitude of moon in radians in its orbit reckoned from A
    σ = s - ξ
    # l (lowercase el) Longitude of moon in its orbit reckoned from its ascending intersection with the equator
    l = σ + 2 * e * np.sin(s - p) + (5. / 4) * e * e * np.sin(2 * (s - p)) + (15. / 4) * m * e * np.sin(s - 2 * h + p)\
        + (11. / 8) * m * m * np.sin(2 * (s - h))

    #
    # Solar Calculations
    #

    # p1 (p-one) Longitude of solar perigee (Schureman [1941, pp. 162])
    # p1 = 281° 13' 15.0" + 6189.03" T + 1.63" T2 + 0.012" T3
    p1 = 4.90822941839 + 0.0300025492114 * T + 7.85398163397e-06 * T2 + 5.3329504922e-08 * T3
    # e1 (e-one) Eccentricity of the Earth's orbit
    e1 = 0.01675104 - 0.00004180 * T - 0.000000126 * T2
    # χ1 (chi-one) right ascension of meridian of place of observations reckoned from the vernal equinox
    χ1 = t + h
    # l1 (lowercase-el(L) one) Longitude of sun in the ecliptic reckoned from the vernal equinox
    l1 = h + 2 * e1 * np.sin(h - p1)
    # cosθ (theta) θ represents the zenith angle of the moon
    cosθ = sinλ * np.sin(I) * np.sin(l) + cosλ * (np.cos(0.5 * I) ** 2 * np.cos(l - χ) + np.sin(0.5 * I) ** 2 *
                                                  np.cos(l + χ))
    # cosφ (phi) φ represents the zenith angle of the run
    cosφ = sinλ * sin(ω) * np.sin(l1) + cosλ * (np.cos(0.5 * ω) ** 2 * np.cos(l1 - χ1) + np.sin(0.5 * ω) ** 2 *
                                                np.cos(l1 + χ1))

    #
    # Distance Calculations
    #

    # (C) Distance parameter, equation 34
    # C**2 = 1/(1 + 0.006738 sinλ ** 2)
    C = np.sqrt(1. / (1 + 0.006738 * sinλ ** 2))
    # (r) Distance from point P to the center of the Earth
    r = C * a + H
    # a' (a prime) Distance parameter, equation 31
    aprime = 1. / (c * (1 - e * e))
    # a1' (a-one prime) Distance parameter, equation 31
    aprime1 = 1. / (c1 * (1 - e1 * e1))
    # (d) Distance between centers of the Earth and the moon
    d = 1. / ((1. / c) + aprime * e * np.cos(s - p) + aprime * e ** 2 * np.cos(2 * (s - p)) + (15. / 8) * aprime * m *
              e * np.cos(s - 2 * h + p) + aprime * m * m * np.cos(2 * (s - h)))
    # (D) Distance between centers of the Earth and the sun
    D = 1. / ((1. / c1) + aprime1 * e1 * np.cos(h - p1))

    # (gm) Vertical component of tidal acceleration due to the moon, equation (1):
    # gm = μMr/d^3 (3 * cos^2(θ) - 1) + 3/2 μMr/d^4 * (5 cos^3 θ - 3 cos(θ))
    gm = (μ * M * r / d ** 3) * (3 * cosθ ** 2 - 1) + (1.5 * (μ * M * r ** 2 / d ** 4) * (5 * cosθ ** 3 - 3 * cosθ))
    # (gs) Vertical component of tidal acceleration due to the sun
    gs = μ * S * r / D ** 3 * (3 * cosφ ** 2 - 1)

    # g0 Total vertical component due to Lunar and Solar forces
    g0 = (gm + gs) * 1e3 * love

    # Returns Lunar, Solar, Total corrections in mGals, as numpy 1d-arrays
    return gm * 1e3 * love, gs * 1e3 * love, g0


def solve_longman_tide_scalar(lat: float, lon: float, alt: float, time: datetime):
    """
    Simple wrapper around solve_longman_tide that allows passing of singular scalar values instead of an array.
    This function simply wraps the scalar values within an ndarray and then extracts the result elements.

    Parameters
    ----------
    lat : float
        Latitude as a scalar float value
    lon : float
        Longitude as a scalar float value
    alt : float
        Altitude as a scalar float value
    time : datetime
        Time as a scalar datetime object

    Returns
    -------
    gm, gs, g0 : float, float, float
        Where gm is the gravitational effect due to the moon
              gs is the gravitational effect due to the sun
              g0 is the total gravitational effect of the sun and moon (gm+gs)

    """
    lat_arr = np.array([lat])
    lon_arr = np.array([lon])
    alt_arr = np.array([alt])
    time_arr = np.array([time])

    gm, gs, g0 = solve_longman_tide(lat_arr, lon_arr, alt_arr, time_arr)
    return gm[0], gs[0], g0[0]


def solve_point_corr(lat: float, lon: float, alt: float, t0=datetime.now(tz=timezone.utc), n=3600, increment='S'):
    """
    Utility function to generate a tide correction DataFrame for a static lat/lon/alt given start time t0,
    an increment, and count (n) of datapoints to generate.
    Default parameters are supplied that will generate a correction series with one second increment over a one hour
    period, with start time being the time of execution.

    Parameters
    ----------
    lat : float
        Latitude in decimal degrees
    lon : float
        Longitude in decimal degrees
    alt : float
        Altitude (height) above sea level in meters
    t0 : DateTime
        Starting date/time
    n : int
        Number of data points to generate
    increment : String
        Increment between data points, uses Pandas offset aliases:
        http://pandas.pydata.org/pandas-docs/stable/timeseries.html#offset-aliases
        Common options:
            H : Hourly Frequency
            T, min : Minutely frequency
            S : Secondly frequency
            L, ms : millisecond frequency
            U, us : microsecond frequency
            N : nanosecond frequency

    Returns
    -------
    df : pd.DataFrame
        DataFrame (index=DateTime, lat, lon, alt, gm, gs, g0) shape: (n, 4)
        Pandas DataFrame indexed by DateTime containing the latitude, longitude, altitude, lunar, solar,
        and total corrections.

    Notes
    -----
    total_corr column renamed to g0 and gm/gs added in v0.3.0

    """
    df = pd.DataFrame(data={'lat': np.repeat(lat, n), 'lon': np.repeat(lon, n), 'alt': np.repeat(alt, n)},
                      index=pd.DatetimeIndex(start=t0, freq=increment, periods=n))

    gm, gs, g0 = solve_longman_tide(df.lat, df.lon, df.alt, df.index)
    df['gm'] = gm
    df['gs'] = gs
    df['g0'] = g0
    return df


def solve_tide_df(df: pd.DataFrame, lat='lat', lon='lon', alt='alt'):
    """
    Solve tidal gravity corrections for a given Pandas DataFrame, returning the source
    DataFrame with an appended 'tide_corr' column.
    Source DataFrame should be indexed by datetime and must contain latitude, longitude,
    and altitude columns.
    Source column names can be specified by passing the applicable column name to the
    lat, lon, alt parameters as required.
    Time is assumed to be the index of the DataFrame

    Returns
    -------
    df : pd.DataFrame
        Copy of source df with added columns 'gm' 'gs' and 'g0': lunar, solar and total corrections respectively

    """
    _lat = df[lat].values
    _lon = df[lon].values
    _alt = df[alt].values
    _time = df.index
    res_df = df.copy(deep=True)  # type: pd.DataFrame

    gm, gs, g0 = solve_longman_tide(_lat, _lon, _alt, _time)
    res_df['gm'] = gm
    res_df['gs'] = gs
    res_df['g0'] = g0
    return res_df

if __name__ == "__main__":
    gdate = datetime(2019, 11, 10, 10, 00, 00)   
    lat = 45
    lon = 105
    alt = 0.0
    gm, gs, g0  = solve_longman_tide_scalar(lat, lon, alt, gdate)
    print(g0)