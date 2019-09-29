/*-----------------------------------------------------------------------------
 *
 * Solar position - original algorithms
 *
 * Author: Martin Paces <martin.paces@eox.at>
 *
 *-----------------------------------------------------------------------------
 * Copyright (C) 2017 EOX IT Services GmbH
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

#ifndef PYSUNPOS_ORIGINAL_H
#define PYSUNPOS_ORIGINAL_H

#include "pymm_aux.h"
#include "sunpos.h"
#include <math.h>

/*---------------------------------------------------------------------------*/
/* recursive array evaluation */

#ifndef RAD2DEG
#define RAD2DEG (180.0/M_PI)
#endif

#ifndef DEG2RAD
#define DEG2RAD (M_PI/180.0)
#endif

static void _sunpos_original(
    ARRAY_DATA arrd_dim,
    ARRAY_DATA arrd_mjd,
    ARRAY_DATA arrd_lat,
    ARRAY_DATA arrd_lon,
    ARRAY_DATA arrd_rad,
    ARRAY_DATA arrd_dtt,
    ARRAY_DATA arrd_pres,
    ARRAY_DATA arrd_temp,
    ARRAY_DATA arrd_out
)
{
    if (arrd_dim.ndim > 0)
    {
        npy_intp i;
        for(i = 0; i < arrd_dim.dim[0]; ++i)
            _sunpos_original(
                _get_arrd_item(&arrd_dim, i),
                _get_arrd_item(&arrd_mjd, i),
                _get_arrd_item(&arrd_lat, i),
                _get_arrd_item(&arrd_lon, i),
                _get_arrd_item(&arrd_rad, i),
                _get_arrd_item(&arrd_dtt, i),
                _get_arrd_item(&arrd_pres, i),
                _get_arrd_item(&arrd_temp, i),
                _get_arrd_item(&arrd_out, i)
            );
    }
    else
    {
        #define P(a,i) ((double*)((char*)(a).data+(i)*(a).stride[0]))
        #define V(a,i) (*P(a,i))

        double rasc, decl, hang, zenith, azimuth;
        double mjd2k, dtt, lat, lon, rad, pres, temp;
        double hours;
        int year, month, day;

        mjd2k = V(arrd_mjd,0); // decimal days since 2000-01-01T00:00:00
        dtt = V(arrd_dtt,0); // seconds (TT to UT offset)
        lat = V(arrd_lat,0) * DEG2RAD; // deg --> rad
        lon = V(arrd_lon,0) * DEG2RAD; // deg --> rad
        rad = V(arrd_rad,0); // km
        pres = V(arrd_pres,0); // atm
        temp = V(arrd_temp,0); // dgC

        mjd2k_to_date(&year, &month, &day, &hours, mjd2k);

        sunpos5original(
            &decl, &rasc, &hang, &azimuth, &zenith,
            year, month, day, hours, lat, lon, rad, dtt, pres, temp
        );

        V(arrd_out,0) = RAD2DEG * decl; // rad --> deg
        V(arrd_out,1) = RAD2DEG * rasc; // rad --> deg
        V(arrd_out,2) = RAD2DEG * hang; // rad --> deg
        V(arrd_out,3) = RAD2DEG * azimuth; // rad --> deg
        V(arrd_out,4) = RAD2DEG * zenith; // rad --> deg

        #undef V
        #undef P
    }
}

#define DOC_SUNPOS_ORIGINAL "\n"\
"   arr_out = sunpos_original(time_mjd2k, lat, lon, rad, dtt, pres, temp)\n\n"\
"     Output:\n"\
"       arr_out - array of the Sun equatorial and horizontal coordinates:\n"\
"                 - declination\n"\
"                 - right ascension\n"\
"                 - hour angle\n"\
"                 - azimuth\n"\
"                 - zenith\n"\
"                 All angles are in deg.\n\n"\
"     Parameters (arrays up to 15 dimensions):\n"\
"       time_mjd2k - array of MJD2000 times.\n"\
"       lat - array of latitudes [deg]\n"\
"       lon - array of longitudes [deg]\n"\
"       rad - array of radii [km] (parallax correction)\n"\
"       dtt - array of offsets to TT [sec] \n"\
"       pres - array of offsets of pressures [atm] (refraction correction) \n"\
"       temp - array of offsets to temperatures [dgC] (refraction coorection)\n"

static PyObject* pysunpos_sunpos_original(PyObject *self, PyObject *args, PyObject *kwdict)
{
    static char *keywords[] = {
        "time_mjd2k", "lat", "lon", "rad", "dtt", "pres", "temp", NULL
    };
    PyObject *obj_mjd = NULL, *arr_mjd = NULL;
    PyObject *obj_lat = NULL, *arr_lat = NULL;
    PyObject *obj_lon = NULL, *arr_lon = NULL;
    PyObject *obj_rad = NULL, *arr_rad = NULL;
    PyObject *obj_dtt = NULL, *arr_dtt = NULL;
    PyObject *obj_pres = NULL, *arr_pres = NULL;
    PyObject *obj_temp = NULL, *arr_temp = NULL;
    PyObject *arr_out = NULL; // return object
    PyObject *retval = NULL; // return object
    PyObject *arr = NULL;
    int idx = 0;

    // parse input arguments
    if (!PyArg_ParseTupleAndKeywords(
        args, kwdict, "OOOOOOO:sunpos", keywords,
        &obj_mjd, &obj_lat, &obj_lon, &obj_rad, &obj_dtt, &obj_pres, &obj_temp
    ))
        goto exit;

    // cast the objects to arrays
    if (NULL == (arr_mjd = _get_as_double_array(obj_mjd, 0, 0, NPY_ALIGNED, keywords[0])))
        goto exit;

    if (NULL == (arr_lat = _get_as_double_array(obj_lat, 0, 0, NPY_ALIGNED, keywords[0])))
        goto exit;

    if (NULL == (arr_lon = _get_as_double_array(obj_lon, 0, 0, NPY_ALIGNED, keywords[0])))
        goto exit;

    if (NULL == (arr_rad = _get_as_double_array(obj_rad, 0, 0, NPY_ALIGNED, keywords[0])))
        goto exit;

    if (NULL == (arr_dtt = _get_as_double_array(obj_dtt, 0, 0, NPY_ALIGNED, keywords[0])))
        goto exit;

    if (NULL == (arr_pres = _get_as_double_array(obj_pres, 0, 0, NPY_ALIGNED, keywords[0])))
        goto exit;

    if (NULL == (arr_temp = _get_as_double_array(obj_temp, 0, 0, NPY_ALIGNED, keywords[0])))
        goto exit;

    // check the array dimensions and allocate the output array
        // get array with the max.dimension
    arr = arr_mjd;
    idx = 0;
    if (PyArray_NDIM(arr) < PyArray_NDIM(arr_lat)) {idx = 1; arr = arr_lat;}
    if (PyArray_NDIM(arr) < PyArray_NDIM(arr_lon)) {idx = 2; arr = arr_lon;}
    if (PyArray_NDIM(arr) < PyArray_NDIM(arr_rad)) {idx = 3; arr = arr_rad;}
    if (PyArray_NDIM(arr) < PyArray_NDIM(arr_dtt)) {idx = 4; arr = arr_dtt;}
    if (PyArray_NDIM(arr) < PyArray_NDIM(arr_pres)) {idx = 5; arr = arr_pres;}
    if (PyArray_NDIM(arr) < PyArray_NDIM(arr_temp)) {idx = 6; arr = arr_temp;}

    // check maximum allowed input array dimension
    if (PyArray_NDIM(arr) > (MAX_OUT_ARRAY_NDIM-1))
    {
        PyErr_Format(PyExc_ValueError, "Array dimension of '%s'"\
            " %d exceeds the allowed maximum value %d!", keywords[idx],
            PyArray_NDIM(arr), (MAX_OUT_ARRAY_NDIM-1));
        goto exit;
    }

    // check array dimensions
    if (PyArray_NDIM(arr_mjd) > 0)
        if (_check_equal_shape(arr_mjd, arr, keywords[0], keywords[idx]))
            goto exit;

    if (PyArray_NDIM(arr_lat) > 0)
        if (_check_equal_shape(arr_lat, arr, keywords[1], keywords[idx]))
            goto exit;

    if (PyArray_NDIM(arr_lon) > 0)
        if (_check_equal_shape(arr_lon, arr, keywords[2], keywords[idx]))
            goto exit;

    if (PyArray_NDIM(arr_rad) > 0)
        if (_check_equal_shape(arr_rad, arr, keywords[3], keywords[idx]))
            goto exit;

    if (PyArray_NDIM(arr_dtt) > 0)
        if (_check_equal_shape(arr_dtt, arr, keywords[4], keywords[idx]))
            goto exit;

    if (PyArray_NDIM(arr_pres) > 0)
        if (_check_equal_shape(arr_pres, arr, keywords[5], keywords[idx]))
            goto exit;

    if (PyArray_NDIM(arr_temp) > 0)
        if (_check_equal_shape(arr_temp, arr, keywords[6], keywords[idx]))
            goto exit;

    // create output array (adding one extra dimension)
    {
        npy_intp *dims_src = PyArray_DIMS(arr);
        npy_intp dims_dst[MAX_OUT_ARRAY_NDIM];
        int i, ndim_src = PyArray_NDIM(arr);

        for (i = 0; i < ndim_src; ++i) dims_dst[i] = dims_src[i];
        dims_dst[ndim_src] = 5;

        if (NULL == (arr_out = PyArray_EMPTY(ndim_src + 1, dims_dst, NPY_DOUBLE, 0)))
            goto exit;
    }

    // evaluate solar positions
    _sunpos_original(
        _array_to_arrd(arr),
        _array_to_arrd(arr_mjd),
        _array_to_arrd(arr_lat),
        _array_to_arrd(arr_lon),
        _array_to_arrd(arr_rad),
        _array_to_arrd(arr_dtt),
        _array_to_arrd(arr_pres),
        _array_to_arrd(arr_temp),
        _array_to_arrd(arr_out)
    );

    // assign return value
    retval = arr_out;

  exit:

    // decrease reference counters to the arrays
    if (arr_mjd) Py_DECREF(arr_mjd);
    if (arr_lat) Py_DECREF(arr_lat);
    if (arr_lon) Py_DECREF(arr_lon);
    if (arr_rad) Py_DECREF(arr_rad);
    if (arr_dtt) Py_DECREF(arr_dtt);
    if (arr_pres) Py_DECREF(arr_pres);
    if (arr_temp) Py_DECREF(arr_temp);
    if (!retval && arr_out) Py_DECREF(arr_out);

    return retval;
}

#endif /* PYSUNPOS_ORIGINAL_H */
