/*-----------------------------------------------------------------------------
 *
 * Geomagnetic Model - C python bindings - coordinate systems definitions
 *
 * Author: Martin Paces <martin.paces@eox.at>
 *
 *-----------------------------------------------------------------------------
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
 *-----------------------------------------------------------------------------
*/

#ifndef PYMM_COORD_H
#define PYMM_COORD_H

typedef enum {
    CT_INVALID = -1,
    CT_GEODETIC_ABOVE_WGS84 = 0,
    CT_GEOCENTRIC_SPHERICAL = 1,
    CT_GEOCENTRIC_CARTESIAN = 2
} COORD_TYPE;

/*
 * Check the coordinate type.
 */
static COORD_TYPE _check_coord_type(int ct, const char *label)
{
    switch (ct)
    {
        case CT_GEODETIC_ABOVE_WGS84:
            return CT_GEODETIC_ABOVE_WGS84;
        case CT_GEOCENTRIC_SPHERICAL:
            return CT_GEOCENTRIC_SPHERICAL;
        case CT_GEOCENTRIC_CARTESIAN:
            return CT_GEOCENTRIC_CARTESIAN;
        default:
            PyErr_Format(PyExc_ValueError, "Invalid coordinate type '%s'!", label);
            return CT_INVALID;
    }
}

#endif  /* PYMM_COORD_H */
