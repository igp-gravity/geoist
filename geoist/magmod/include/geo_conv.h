/**
 * @file geo_conv.h
 * @author Martin Paces <martin.paces@eox.at>
 * @brief Geo-coordinates conversions.
 *
 * This file contains various geo-coordinates conversion utilities.
 *
 * Copyright (C) 2014 EOX IT Services GmbH
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies of this Software or works derived from this Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

#ifndef GEO_CONV_H
#define GEO_CONV_H 1

#include <math.h>
#include "math_aux.h"

#define DG2RAD (M_PI/180.0)
#define RAD2DG (180.0*M_1_PI)

#define WGS84_A 6378.137
//#define WGS84_B 6356.752314245
#define WGS84_B 6356.7523142
#define WGS84_INVFLAT 298.257223563
#define WGS84_FLAT (1.0/WGS84_INVFLAT)
#define WGS84_EPS2 (1.0 - (WGS84_B*WGS84_B)/(WGS84_A*WGS84_A))
#define WGS84_EPS sqrt(WGS84_EPS2)
#define WGS84_RADIUS 6371.2
#ifndef M_PI
#define M_PI 3.14159265358979323846  /* pi */
#endif

#ifndef M_1_PI
#define M_1_PI 0.31830988618379067154  /* 1/pi */
#endif

/**
 * @brief Low level conversion from geodetic to geocentric coordinates.
 *
 * Required ellipsoid parameters:
 *  elp_a - semi-major axis
 *  elp_e2 - squared first eccentricity
 */

static void _geodetic2geocentric(
    double *rad_xy, double *z,
    double lat_rad, double elv,
    double elp_a, double elp_e2)
{
    const double sin_lat = sin(lat_rad);
    const double cos_lat = cos(lat_rad);

    /* radius of curvature of the ellipsoid */
    const double rc = elp_a / sqrt(1.0 - elp_e2*sin_lat*sin_lat);

    /* radius on the equator plane */
    *rad_xy = (rc + elv)*cos_lat;

    /* z coordinate */
    *z = (rc*(1.0 - elp_e2) + elv)*sin_lat;

}

/**
 * @brief Low level conversion from geocentric to geodetic geocentric coordinates.
 *
 * Uses Ferrari's solution.
 *
 * Required ellipsoid parameters:
 *  elp_a - semi-major axis
 *  elp_e2 - squared first eccentricity
 */

static void _geocentric2geodetic(
    double *lat_rad, double *elv,
    double rad_xy, double z,
    double elp_a, double elp_e2)
{
    const double p = rad_xy;
    /* Ferrari's solution */
    const double pa = p/elp_a;
    const double za = z/elp_a;
    const double pa2 = pa*pa;
    const double za2 = za*za;
    const double ee4 = elp_e2*elp_e2;
    const double rkp0 = (1.0 - elp_e2);
    const double zt = rkp0*za2;
    const double rh = (pa2 + zt - ee4)/6.0;
    const double ss = 0.25*zt*ee4*pa2;
    const double rh3 = rh*rh*rh;
    const double tmp = rh3 + ss + sqrt(ss*(ss+2.0*rh3));
    const double tt = copysign(pow(fabs(tmp), 1.0/3.0), tmp);
    const double uu = rh + tt + rh*rh/tt;
    const double vv = sqrt(uu*uu + ee4*zt);
    const double ww = 0.5*elp_e2*(uu + vv - zt)/vv;
    const double kp = 1.0 + elp_e2*(sqrt(uu + vv + ww*ww) + ww)/(uu + vv);
    const double zkp = z * kp;

    *lat_rad = atan2(zkp, p);
    *elv = norm2d(p, zkp)*(1.0/kp - rkp0)/elp_e2;
}


/**
 * @brief Convert geodetic to geocentric Cartesian coordinates.
 *
 * Covert geodetic coordinates (latitude, longitude, elevation(above ellipsoid)
 * to geocentric Cartesian (ECEF) coordinates.
 *
 * The input geodetic coordinates shall be in degrees. The units of the output
 * Cartesian coordinates are the same as the one of the ellipsoid semi-major axis.
 *
 * Required ellipsoid parameters:
 *  elp_a - semi-major axis
 *  elp_e2 - squared first eccentricity
 */

static void geodetic2geocentric_cart(
    double *x, double *y, double *z,
    double lat, double lon, double elv,
    double elp_a, double elp_e2)
{
    const double lat_rad = DG2RAD*lat;
    const double lon_rad = DG2RAD*lon;
    double rad_xy;

    _geodetic2geocentric(&rad_xy, z, lat_rad, elv, elp_a, elp_e2);

    *x = rad_xy * cos(lon_rad);
    *y = rad_xy * sin(lon_rad);
}


/**
 * @brief Convert geocentric coordinates to geocentric spherical coordinates
 *
 * Covert geodetic coordinates (latitude, longitude, elevation(above ellipsoid)
 * to geocentric spherical coordinates (radius, elevation/theta, azimuth/phi).
 *
 * The input geodetic coordinates shall be in degrees. The units of the output
 * Cartesian coordinates as well as the unit of the radial coordinate are
 * the same as the one of the ellipsoid semi-major axis. The output spherical
 * coordinates are in radians.
 *
 * Required ellipsoid parameters:
 *  elp_a - semi-major axis
 *  elp_e2 - squared first eccentricity
 */

static void geodetic2geocentric_sph(
    double *r, double *th, double *ph,
    double lat, double lon, double elv,
    double elp_a, double elp_e2)
{
    const double lat_rad = DG2RAD*lat;
    const double lon_rad = DG2RAD*lon;
    double rad_xy, z;

    _geodetic2geocentric(&rad_xy, &z, lat_rad, elv, elp_a, elp_e2);

    *r = norm2d(rad_xy, z);
    *th = atan2(z, rad_xy);
    *ph = lon_rad;
}


/**
 * @brief Convert geodetic coordinates to geocentric Cartesian and spherical ones.
 *
 * Covert geodetic coordinates (latitude, longitude, elevation(above ellipsoid))
 * to geocentric coordinates both Cartesian (x, y, z) and spherical (radius,
 * elevation/theta, azimuth/phi). The input geodetic coordinates shall be in
 * degrees. The output spherical coordinates are in radians.
 * The input geodetic coordinates shall be in degrees. The unit of the output
 * the radial coordinate is  the same as the one of the ellipsoid semi-major
 * axis. The output spherical coordinates are in radians.
 *
 * Use when both, spherical and Cartesian, coordinates are needed.
 * Note that this function is more efficient that two calls to
 *  geodetic2geocentric_cart
 *  geodetic2geocentric_sph
 *
 * Required ellipsoid parameters:
 *  elp_a - semi-major axis
 *  elp_e2 - squared first eccentricity
 */

static void geodetic2geocentric(
    double *x, double *y, double *z,
    double *r, double *th, double *ph,
    double lat, double lon, double elv,
    double elp_a, double elp_e2)
{
    const double lat_rad = DG2RAD*lat;
    const double lon_rad = DG2RAD*lon;
    double rad_xy;

    _geodetic2geocentric(&rad_xy, z, lat_rad, elv, elp_a, elp_e2);

    *x = rad_xy * cos(lon_rad);
    *y = rad_xy * sin(lon_rad);

    *r = norm2d(rad_xy, *z);
    *th = atan2(*z, rad_xy);
    *ph = lon_rad;
}

/**
 * @brief Convert spherical to Cartesian coordinates
 *
 * Covert spherical coordinates(radius, elevation/theta, azimuth/phi)
 * to Cartesian coordinates (x, y, z).
 *
 * The spherical coordinates shall be in radians.
 */

static void sph2cart(
    double *x, double *y, double *z,
    double r, double th, double ph)
{
    const double sin_th = sin(th);
    const double cos_th = cos(th);
    const double sin_ph = sin(ph);
    const double cos_ph = cos(ph);

    /* radius on the azimuth plane (z=0)*/
    const double ra = r*cos_th;

    *x = ra*cos_ph;
    *y = ra*sin_ph;
    *z = r*sin_th;
}

/**
 * @brief Convert Cartesian to spherical coordinates
 *
 * Covert Cartesian coordinates (x, y, z) to spherical
 * coordinates (radius, elevation/theta, azimuth/phi)
 *
 * The spherical coordinates are produced in radians.
 */
static void cart2sph(
    double *r, double *th, double *ph,
    double x, double y, double z)
{
    *r = norm3d(x, y, z);
    *th = asin(z/(*r));
    *ph = atan2(y, x);
}


/**
 * @brief Convert Cartesian geocentric to geodetic coordinates
 *
 * Covert Cartesian (ECEF) coordinates (x, y, z) to geodetic coordinates
 * (latitude, longitude, elevation(above ellipsoid).
 *
 * The geodetic coordinates are produced in degrees. The unit of the elevation
 * is the same as the unit of the Cartesian coordinates and of the ellipsoid
 * semi-major axis.
 */

static void geocentric_cart2geodetic(
    double *lat, double *lon, double *elv,
    double x, double y, double z,
    double elp_a, double elp_e2)
{
    double lat_rad;

    _geocentric2geodetic(&lat_rad, elv, norm2d(x, y), z, elp_a, elp_e2);

    *lat = RAD2DG*lat_rad;
    *lon = RAD2DG*atan2(y, x);
}

/**
 * @brief Convert spherical geocentric to geodetic coordinates
 *
 * Convert geocentric spherical coordinates (radius, elevation/theta, azimuth/phi)
 * to geodetic coordinates (latitude, longitude, elevation(above ellipsoid).
 *
 * The spherical coordinates are expect to be in radians.
 * The geodetic coordinates are produced in degrees. The unit of the elevation
 * is the same as the unit of the radial coordinate and of the ellipsoid
 * semi-major axis.
 */

static void geocentric_sph2geodetic(
    double *lat, double *lon, double *elv,
    double r, double th, double ph,
    double elp_a, double elp_e2)
{
    double lat_rad;

    _geocentric2geodetic(&lat_rad, elv, r*cos(th), r*sin(th), elp_a, elp_e2);

    *lat = RAD2DG*lat_rad;
    *lon = RAD2DG*ph;
}

#endif  /*GEO_CONV_H*/
