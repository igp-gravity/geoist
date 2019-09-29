/*-----------------------------------------------------------------------------
 *
 * Geomagnetic Model - C python bindings - vector rotation - cart2sph
 *  (i.e., vector coordinate system transformation)
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

#ifndef PYMM_VROT_CART2SPH_H
#define PYMM_VROT_CART2SPH_H

#include <math.h>
#include "math_aux.h"
#include "geo_conv.h"
#include "pymm_aux.h"
#include "pymm_vrot_common.h"

/* recursive vector rotation */

static void _vrot_cart2sph(ARRAY_DATA ad_i, ARRAY_DATA ad_lat,
                           ARRAY_DATA ad_lon, ARRAY_DATA ad_o)
{
    if (ad_i.ndim > 1)
    {
        npy_intp i;
        for(i = 0; i < ad_i.dim[0]; ++i)
            _vrot_cart2sph(
                _get_arrd_item(&ad_i, i),
                _get_arrd_item(&ad_lat, i),
                _get_arrd_item(&ad_lon, i),
                _get_arrd_item(&ad_o, i)
            );
    }
    else
    {
        #define SV(a) (*((double*)((a).data)))
        #define P(a,i) ((double*)((char*)(a).data+(i)*(a).stride[0])) //((double*)((a).data+(i)*(a).stride[0]))
        #define V(a,i) (*P(a,i))

        const double lat = -DG2RAD*SV(ad_lat);
        const double lon = -DG2RAD*SV(ad_lon);
        double tmp;

        // rotate around the azimuth axis
        rot2d(&tmp, P(ad_o,1), V(ad_i,0), V(ad_i,1), sin(lon), cos(lon));

        // rotate around the elevation axis
        rot2d(P(ad_o,2), P(ad_o,0), tmp, V(ad_i,2), sin(lat), cos(lat));

        #undef V
        #undef P
    }
}

/* python function definition */

#define DOC_VROT_CART2SPH "\n"\
"   arr_out = vrot_cart2sph(arr_in, arr_lat, arr_lon)\n"\
"\n"\
"     Rotate vectors from the Cartesian (XYZ) to \n"\
"     the spherical (NEC) coordinate frame for the given latitude\n"\
"     and longitude in degrees.\n"\
"\n"\
"     The inputs are:\n"\
"         arr_in - array of the input vectors\n"\
"         arr_lat - array of latitudes.\n"\
"         arr_lon - array of longitudes.\n"\
"     Scalar lat/lon values are also accepted for a single vector rotation.\n"

static PyObject* vrot_cart2sph(PyObject *self, PyObject *args, PyObject *kwdict)
{
    static char *keywords[] = {"arr_in", "arr_lat", "arr_lon", NULL};
    PyObject *obj_in = NULL; // input object
    PyObject *obj_lat = NULL; // input object
    PyObject *obj_lon = NULL; // input object
    PyObject *arr_in = NULL; // input array
    PyObject *arr_lat = NULL; // input array
    PyObject *arr_lon = NULL; // input array
    PyObject *arr_out = NULL; // output array
    PyObject *retval = NULL;

    // parse input arguments
    if (!PyArg_ParseTupleAndKeywords(
        args, kwdict, "OOO|:vrot_cart2sph", keywords,
        &obj_in, &obj_lat, &obj_lon
    ))
        goto exit;

    // cast the objects to arrays
    if (NULL == (arr_in=_get_as_double_array(obj_in, 1, 0, NPY_ALIGNED, keywords[0])))
        goto exit;

    if (NULL == (arr_lat=_get_as_double_array(obj_lat, 0, 0, NPY_ALIGNED, keywords[1])))
        goto exit;

    if (NULL == (arr_lon=_get_as_double_array(obj_lon, 0, 0, NPY_ALIGNED, keywords[2])))
        goto exit;

    // check maximum allowed input array dimension
    if (PyArray_NDIM(arr_in) > MAX_OUT_ARRAY_NDIM)
    {
        PyErr_Format(PyExc_ValueError, "The input array dimension of '%s'"\
            " %d exceeds the allowed maximum value %d!", keywords[0],
            PyArray_NDIM(arr_in), MAX_OUT_ARRAY_NDIM);
        goto exit;
    }

    // check the dimensions
    if (_check_array_dim_eq(arr_in, -1, 3, keywords[0]))
        goto exit;

    if (_vrot_arr_check(arr_in, arr_lat, keywords[0], keywords[1]))
        goto exit;

    if (_vrot_arr_check(arr_in, arr_lon, keywords[0], keywords[2]))
        goto exit;

    // create the output array
    if (NULL == (arr_out = _get_new_double_array(PyArray_NDIM(arr_in), PyArray_DIMS(arr_in), 3)))
        goto exit;

    // rotate the vector(s)
    _vrot_cart2sph(_array_to_arrd(arr_in), _array_to_arrd(arr_lat),
                   _array_to_arrd(arr_lon), _array_to_arrd(arr_out));

    retval = arr_out;

  exit:

    // decrease reference counters of the arrays
    if (arr_in){Py_DECREF(arr_in);}
    if (arr_lat){Py_DECREF(arr_lat);}
    if (arr_lon){Py_DECREF(arr_lon);}
    if (!retval && arr_out){Py_DECREF(arr_out);}

    return retval;
}

#endif  /* PYMM_VROT_CART2SPH_H */
