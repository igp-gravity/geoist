/*-----------------------------------------------------------------------------
 *
 * Geomagnetic Model - C python bindings - vector rotation - sph2geod
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

#ifndef PYMM_VROT_SPH2GEOD_H
#define PYMM_VROT_SPH2GEOD_H

#include "math_aux.h"
#include "geo_conv.h"
#include "pymm_aux.h"
#include "pymm_vrot_common.h"

/* recursive vector rotation */

static void _vrot_sph2geod(ARRAY_DATA ad_i, ARRAY_DATA ad_dlat,
                           ARRAY_DATA ad_o)
{
    if (ad_i.ndim > 1)
    {
        npy_intp i;
        for(i = 0; i < ad_i.dim[0]; ++i)
            _vrot_sph2geod(
                _get_arrd_item(&ad_i, i),
                _get_arrd_item(&ad_dlat, i),
                _get_arrd_item(&ad_o, i)
            );
    }
    else
    {
        #define SV(a) (*((double*)((a).data)))
        #define P(a,i) ((double*)((char*)(a).data+(i)*(a).stride[0])) //((double*)((a).data+(i)*(a).stride[0]))
        #define V(a,i) (*P(a,i))

        const double dlat = -DG2RAD*SV(ad_dlat);

        // rotate around the elevation axis
        rot2d(P(ad_o,2), P(ad_o,0), V(ad_i,2), V(ad_i,0), sin(dlat), cos(dlat));
        V(ad_o,1) = V(ad_i,1);

        #undef V
        #undef P
    }
}

/* python function definition */

#define DOC_VROT_SPH2GEOD "\n"\
"   arr_out = vrot_sph2geod(arr_in, arr_dlat)\n"\
"\n"\
"     Rotate vectors from the geocentric spherical (NEC) to the geodetic\n"\
"     (NEC) coordinate frame for a given difference of latitudes\n"\
"     in degrees.\n"\
"     This function can be also used for the inverse rotation from the geodetic\n"\
"     (NEC) to the geocentric spherical (NEC) coordinate frame by setting\n"\
"     the negative difference of latitudes.\n"\
"\n"\
"     The inputs are:\n"\
"         arr_in - array of the input vectors\n"\
"         arr_dlat - array of differences of the latitudes.\n"\
"     A scalar value is also accepted for a single vector rotation.\n"

static PyObject* vrot_sph2geod(PyObject *self, PyObject *args, PyObject *kwdict)
{
    static char *keywords[] = {"arr_in", "arr_dlat", NULL};
    PyObject *obj_in = NULL; // input object
    PyObject *obj_dlat = NULL; // input object
    PyObject *arr_in = NULL; // input array
    PyObject *arr_dlat = NULL; // input array
    PyObject *arr_out = NULL; // output array
    PyObject *retval = NULL;

    // parse input arguments
    if (!PyArg_ParseTupleAndKeywords(
        args, kwdict, "OO|:vrot_sph2geod", keywords, &obj_in, &obj_dlat
    ))
        goto exit;

    // cast the objects to arrays
    if (NULL == (arr_in=_get_as_double_array(obj_in, 1, 0, NPY_ALIGNED, keywords[0])))
        goto exit;

    if (NULL == (arr_dlat=_get_as_double_array(obj_dlat, 0, 0, NPY_ALIGNED, keywords[1])))
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

    if (_vrot_arr_check(arr_in, arr_dlat, keywords[0], keywords[1]))
        goto exit;

    // create the output array
    if (NULL == (arr_out = _get_new_double_array(PyArray_NDIM(arr_in), PyArray_DIMS(arr_in), 3)))
        goto exit;

    // rotate the vector(s)
    _vrot_sph2geod(_array_to_arrd(arr_in), _array_to_arrd(arr_dlat), _array_to_arrd(arr_out));

    retval = arr_out;

  exit:

    // decrease reference counters to the arrays
    if (arr_in){Py_DECREF(arr_in);}
    if (arr_dlat){Py_DECREF(arr_dlat);}
    if (!retval && arr_out){Py_DECREF(arr_out);}

    return retval;
}

#endif  /* PYMM_VROT_SPH2GEOD_H */
