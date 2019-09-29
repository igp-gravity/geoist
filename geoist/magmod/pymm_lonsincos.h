/*-----------------------------------------------------------------------------
 *
 * Geomagnetic Model - C python bindings
 * - longitude sin/cos spherical terms' evaluation
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

#ifndef PYMM_LONSINCOS_H
#define PYMM_LONSINCOS_H

#include <math.h>
#include "shc.h"
#include "pymm_aux.h"

/* python function definition */

#define DOC_LONSINCOS "\n"\
"  lonsin, loncos = lonsincos(longitude, degree, fast_alg=True)\n"\
"\n"\
"     For given 'longitude' (geocentric spherical) and 'degree', evaluate\n"\
"     the following sine and cosine series: \n"\
"        cos(i*longitude) for i in range(0, degree+1)\n"\
"        sin(i*longitude) for i in range(0, degree+1)\n"\
"     The longitude has to be entered in dg..\n"\
"     The 'fast_alg' boolean options forces the subroutine to use a faster\n"\
"     but slightly less precise evaluation algorithm.\n"

static PyObject* lonsincos(PyObject *self, PyObject *args, PyObject *kwdict)
{
    static char *keywords[] = {"longitude", "degree", "fast_alg", NULL};

    int degree;
    int fast_alg = 1;
    double lon_dg;
    PyObject *arr_lonsin = NULL; // lonsin array
    PyObject *arr_loncos = NULL; // loncos array
    PyObject *retval = NULL; // output tuple

    // parse input arguments
    if (!PyArg_ParseTupleAndKeywords(
        args, kwdict, "di|i:lonsincos", keywords, &lon_dg, &degree, &fast_alg
    ))
        goto exit;

    if (degree < 0)
    {
        PyErr_Format(PyExc_ValueError, "%s < 0", keywords[1]);
        goto exit;
    }

    // create the output arrays
    if (NULL == (arr_lonsin = _get_new_double_array(1, NULL, degree+1)))
        goto exit;

    if (NULL == (arr_loncos = _get_new_double_array(1, NULL, degree+1)))
        goto exit;

    if (NULL == (retval = Py_BuildValue("NN", arr_lonsin, arr_loncos)))
        goto exit;

    // evaluate series
    if (fast_alg)
        shc_azmsincos(PyArray_DATA(arr_lonsin), PyArray_DATA(arr_loncos), degree, DG2RAD*lon_dg);
    else
        shc_azmsincos_ref(PyArray_DATA(arr_lonsin), PyArray_DATA(arr_loncos), degree, DG2RAD*lon_dg);

  exit:

    // decrease reference counters to the arrays
    if (!retval && arr_lonsin){Py_DECREF(arr_lonsin);}
    if (!retval && arr_loncos){Py_DECREF(arr_loncos);}

    return retval;
}

#endif  /* PYMM_LONSINCOS_H */
