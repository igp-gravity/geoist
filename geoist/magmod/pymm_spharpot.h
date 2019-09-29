/*-----------------------------------------------------------------------------
 *
 * Geomagnetic Model - C python bindings - spherical-harmonics
 * - evaluation of the potential gradient in the spherical coordinate system
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

#ifndef PYMM_SPHARPOT_H
#define PYMM_SPHARPOT_H

#include "shc.h"
#include "pymm_aux.h"
#include "pymm_cconv.h"

/* python function definition */
#define DOC_SPHARPOT "\n"\
"  v_grad = spharpot(radius, degree, coef_g, coef_h, leg_p, rrp, lonsin, loncos)\n"\
"\n"\
"     Spherical harmonic evaluation of the value of the potential\n"\
"     (scalar) field in the (geocentric) spherical coordinates (latitude,\n"\
"     longitude, radius).\n"\
"     The input parameters are:\n"\
"       radius - radius, i.e., distance from the earth centre, at the\n"\
"                evaluated location.\n"\
"       degree - degree of the spherical harmonic model.\n"\
"       coef_g - vector of spherical harmonic model coefficients.\n"\
"       coef_h - vector of spherical harmonic model coefficients.\n"\
"       leg_p - vector the Legendre polynomials.\n"\
"       rrp - vector the relative radius powers.\n"\
"       lonsin - vector of the longitude cosines.\n"\
"       lonsin - vector of the longitude sines.\n"

static PyObject* spharpot(PyObject *self, PyObject *args, PyObject *kwdict)
{
    static char *keywords[] = {"radius", "degree", "coef_g", "coef_h",
                                "leg_p", "rpp", "lonsin", "loncos", NULL};

    int degree, nterm;
    double rad;

    PyObject *obj_cg = NULL; // coef_g object
    PyObject *obj_ch = NULL; // coef_h object
    PyObject *obj_lp = NULL; // P object
    PyObject *obj_rrp = NULL; // rel.rad.pow. object
    PyObject *obj_lsin = NULL; // lonsin object
    PyObject *obj_lcos = NULL; // loncos object

    PyObject *arr_out = NULL; // output array
    PyObject *arr_cg = NULL; // coef_g array
    PyObject *arr_ch = NULL; // coef_h array
    PyObject *arr_lp = NULL; // P array
    PyObject *arr_rrp = NULL; // rel.rad.pow. array
    PyObject *arr_lsin = NULL; // lonsin array
    PyObject *arr_lcos = NULL; // loncos array

    // parse input arguments
    if (!PyArg_ParseTupleAndKeywords(
        args, kwdict, "diOOOOOO|:spharpot", keywords,
        &rad, &degree, &obj_cg, &obj_ch, &obj_lp, &obj_rrp, &obj_lsin, &obj_lcos
    ))
        goto exit;

    if (degree < 0)
    {
        PyErr_Format(PyExc_ValueError, "%s < 0", keywords[1]);
        goto exit;
    }

    nterm = ((degree+1)*(degree+2))/2;

    // cast the objects to arrays
    if (NULL == (arr_cg=_get_as_double_array(obj_cg, 1, 1, NPY_IN_ARRAY, keywords[2])))
        goto exit;

    if (NULL == (arr_ch=_get_as_double_array(obj_ch, 1, 1, NPY_IN_ARRAY, keywords[3])))
        goto exit;

    if (NULL == (arr_lp=_get_as_double_array(obj_lp, 1, 1, NPY_IN_ARRAY, keywords[4])))
        goto exit;

    if (NULL == (arr_rrp=_get_as_double_array(obj_rrp, 1, 1, NPY_IN_ARRAY, keywords[5])))
        goto exit;

    if (NULL == (arr_lsin=_get_as_double_array(obj_lsin, 1, 1, NPY_IN_ARRAY, keywords[6])))
        goto exit;

    if (NULL == (arr_lcos=_get_as_double_array(obj_lcos, 1, 1, NPY_IN_ARRAY, keywords[7])))
        goto exit;

    // check the arrays' dimensions
    if (_check_array_dim_le(arr_cg, 0, nterm, keywords[2]))
        goto exit;

    if (_check_array_dim_le(arr_ch, 0, nterm, keywords[3]))
        goto exit;

    if (_check_array_dim_le(arr_lp, 0, nterm, keywords[4]))
        goto exit;

    if (_check_array_dim_le(arr_rrp, 0, degree+1, keywords[5]))
        goto exit;

    if (_check_array_dim_le(arr_lsin, 0, degree+1, keywords[6]))
        goto exit;

    if (_check_array_dim_le(arr_lcos, 0, degree+1, keywords[7]))
        goto exit;

    // allocate the output numpy scalar (zero-dimensional array)
    if (NULL == (arr_out = PyArray_EMPTY(0, NULL, NPY_DOUBLE, 0)))
        goto exit;

    // the evaluation
    {
        double *out = (double*)PyArray_DATA(arr_out);
        shc_eval(
            out, NULL, NULL, NULL, degree, 0x1,
            0.0, rad, PyArray_DATA(arr_cg),
            PyArray_DATA(arr_ch), PyArray_DATA(arr_lp),
            PyArray_DATA(arr_lp), PyArray_DATA(arr_rrp),
            PyArray_DATA(arr_lsin), PyArray_DATA(arr_lcos),
            0 // this flag has no impact on the evaluated potential
        );
    }

  exit:

    // decrease reference counters to the arrays
    if (arr_cg){Py_DECREF(arr_cg);}
    if (arr_ch){Py_DECREF(arr_ch);}
    if (arr_lp){Py_DECREF(arr_lp);}
    if (arr_rrp){Py_DECREF(arr_rrp);}
    if (arr_lsin){Py_DECREF(arr_lsin);}
    if (arr_lcos){Py_DECREF(arr_lcos);}

    return arr_out;
 }

#endif  /* PYMM_SPHARPOT_H */
