/*-----------------------------------------------------------------------------
 *
 * Geomagnetic Model - C python bindings
 * - associative Legendre functions evaluation
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

#ifndef PYMM_LEGENDRE_H
#define PYMM_LEGENDRE_H

#include "shc.h"
#include "pymm_aux.h"
#include "pymm_cconv.h"

/* python function definition */

#define DOC_LEGENDRE "\n"\
"   p, dp = legendre(latitude, degree)\n"\
"\n"\
"     For given the 'latitude' in degrees and the model's 'degree' evaluate\n"\
"     the Schmidt semi-normalised associated Legendre functions and its\n"\
"     and its derivatives: \n"\
"         P_n^m(sin(latitude))  and  dP_n^m(sin(latitude))\n"\
"     where n = 0..degree and m = 0..n \n"\
"\n"\
"     The input parameters are:\n"\
"       latitude - spherical latitude in degrees.\n"\
"       degree - degree of the spherical harmonic model.\n"\
"\n"


static PyObject* legendre(PyObject *self, PyObject *args, PyObject *kwdict)
{
    static char *keywords[] = {"latitude", "degree", NULL};

    int degree, nterm;
    double lat_sph;
    PyObject *arr_p = NULL; // P array
    PyObject *arr_dp = NULL; // dP array
    PyObject *retval = NULL; // output tuple

    // parse input arguments
    if (!PyArg_ParseTupleAndKeywords(
        args, kwdict, "di:legendre", keywords, &lat_sph, &degree
    ))
        goto exit;

    if (degree < 0)
    {
        PyErr_Format(PyExc_ValueError, "%s < 0", keywords[1]);
        goto exit;
    }

    if ((lat_sph < -90.0) || (lat_sph > 90.0))
    {
        if (lat_sph < 0.0) {
            PyErr_Format(PyExc_ValueError, "%s < -90", keywords[0]);
        } else {
            PyErr_Format(PyExc_ValueError, "%s > 90", keywords[0]);
        }
        goto exit;
    }

    nterm = ((degree+1)*(degree+2))/2;

    // create the output arrays
    if (NULL == (arr_p = _get_new_double_array(1, NULL, nterm)))
        goto exit;

    if (NULL == (arr_dp = _get_new_double_array(1, NULL, nterm)))
        goto exit;

    {
        // allocate and fill the pre-calculated square-roots
        double *psqrt = NULL;
        if (NULL == (psqrt = shc_presqrt(degree)))
        {
            PyErr_Format(PyExc_ValueError, "Memory allocation error!");
            goto exit;
        }

        shc_legendre(PyArray_DATA(arr_p), PyArray_DATA(arr_dp), degree, DG2RAD*lat_sph, psqrt);

        // free the square root array
        free(psqrt);
    }

    if (NULL == (retval = Py_BuildValue("NN", arr_p, arr_dp)))
        goto exit;

  exit:

    // decrease reference counters to the arrays
    if (!retval && arr_p){Py_DECREF(arr_p);}
    if (!retval && arr_dp){Py_DECREF(arr_dp);}

    return retval;
}

#endif  /* PYMM_LEGENDRE_H */

