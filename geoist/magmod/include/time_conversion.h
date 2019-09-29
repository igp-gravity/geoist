/**
 * @file time_conversion.h
 * @author Martin Paces <martin.paces@eox.at>
 * @brief Time Conversion Subroutines
 *
 * Time conversion subroutines.
 *
 * Copyright (C) 2018 EOX IT Services GmbH
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
#ifndef TIME_CONVERSION_H
#define TIME_CONVERSION_H

#include <math.h>

/**
 * @brief Convert MJD2000 to year, month, day and day fraction
 *
 * Convert MJD2000 to year, month, day and number of decimal hours
 *
 * ref: https://en.wikipedia.org/wiki/Julian_day#Julian_or_Gregorian_calendar_from_Julian_day_number
 * Gregorian date formula applied since 1582-10-15
 * Julian date formula applied until 1582-10-04
 */

static void mjd2k_to_date(int *year, int *month, int *day, double *hours, double mjd2k) {
    const int day2k = (int)floor(mjd2k);
    const int d__ = day2k + 2451545;
    const int f0_ = d__ + 1401 ;
    const int f__ = f0_ + (d__ > 2299160 ? (((4*d__ + 274277)/146097)*3)/4 - 38: 0);
    const int e__ = 4*f__ + 3;
    const int h__ = 5*((e__ % 1461)/4) + 2;
    *day = (h__%153)/5 + 1;
    *month = (h__/153 + 2)%12 + 1;
    *year = e__/1461 - 4716 + (14 - *month)/12;
    *hours = 24.0 * (mjd2k - day2k);
}


/**
 * @brief Convert MJD2000 to a whole year.
 *
 * ref: https://en.wikipedia.org/wiki/Julian_day#Julian_or_Gregorian_calendar_from_Julian_day_number
 * Gregorian date formula applied since 1582-10-15
 * Julian date formula applied until 1582-10-04
 */

static int mjd2k_to_year(double mjd2k) {
    const int day2k = (int)floor(mjd2k);
    const int d__ = day2k + 2451545;
    const int f0_ = d__ + 1401 ;
    const int f__ = f0_ + (d__ > 2299160 ? (((4*d__ + 274277)/146097)*3)/4 - 38: 0);
    const int e__ = 4*f__ + 3;
    const int h__ = 5*((e__ % 1461)/4) + 2;
    const int month = (h__/153 + 2)%12 + 1;
    return e__/1461 - 4716 + (14 - month)/12;
}


/**
 * @brief Get MJD2000 of the year start.
 *
 * The function returns number of days between the start of the given year
 * and start of the year 2000.
 *
 * ref: https://en.wikipedia.org/wiki/Julian_day#Converting_Julian_or_Gregorian_calendar_date_to_Julian_day_number
 * ref: http://www.cs.utsa.edu/~cs1063/projects/Spring2011/Project1/jdn-explanation.html
 * Gregorian date formula applied since 1583
 * Julian date formula applied until 1582
 */

static int mjd2k_year_start(int year) {
    const int y__ = year + 4799;
    return (365*y__ + y__/4 - 2483321 + (year > 1582)*(y__/400 - y__/100 + 38));
}


/**
 * @brief Return 1 for a leap year or 0 otherwise.
 *
 * Gregorian date formula applied since 1583
 * Julian date formula applied until 1582
 */
static int is_leap_year(int year) {
    return ((year%4 == 0)&&((year%100 != 0)||(year%400 == 0)||(year <= 1582)));
}

/**
 * @brief Get number of days in the given year.
 *
 * This function returns 356 for a leap year or 366 otherwise.
 *
 * Gregorian date formula applied since 1583
 * Julian date formula applied until 1582
 */
static int days_per_year(int year) {
    return 365 + is_leap_year(year);
}

/**
 * @brief Get fraction of a year passed in one day.
 *
 * This function returns 1/366 for a leap year or 1/355 otherwise.
 *
 * Gregorian date formula applied since 1583
 * Julian date formula applied until 1582
 */
static double years_per_day(int year) {
    return is_leap_year(year) ? 1./366. : 1./365.;
}

/**
 * @brief Convert MJD2000 to decimal year.
 *
 * Gregorian date formula applied since 1582-10-15
 * Julian date formula applied until 1582-10-04
 */

static double mjd2k_to_decimal_year(double mjd2k) {
    const int year = mjd2k_to_year(mjd2k);
    return year + (mjd2k - mjd2k_year_start(year))*years_per_day(year);
}

/**
 * @brief Convert MJD2000 to year fraction.
 *
 * Equivalent to:
 *  decimal_year = mjd2k_to_decimal_year(mjd2k);
 *  year_fraction = decimal_year - floor(decimal_year);
 */

static double mjd2k_to_year_fraction(double mjd2k) {
    const int year = mjd2k_to_year(mjd2k);
    return (mjd2k - mjd2k_year_start(year))*years_per_day(year);
}

/**
 * @brief Convert decimal year to MJD2000.
 *
 * Gregorian date formula applied since 1582-10-15
 * Julian date formula applied until 1582-10-04
 */

static double decimal_year_to_mjd2k(double decimal_year) {
    if (isfinite(decimal_year)) {
        const double whole_year = floor(decimal_year);
        const int year = (int)whole_year;
        return mjd2k_year_start(year) + (decimal_year - whole_year)*days_per_year(year);
    } else {
        return decimal_year;
    }
}

#endif /* TIME_CONVERSION */
