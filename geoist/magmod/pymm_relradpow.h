/*-----------------------------------------------------------------------------
 *
 * Geomagnetic Model - C python bindings
 * - evaluation of the relative radius power series
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

#ifndef PYMM_RELRADPOW_H
#define PYMM_RELRADPOW_H

#include "shc.h"
#include "pymm_aux.h"

/* Earth radius in km */
#define RADIUS  6371.2

#ifndef STR
#define SRINGIFY(x) #x
#define STR(x) SRINGIFY(x)
#endif

/* python function definition */

#define DOC_RELRADPOW "\n"\
"   rrp = relradpow(radius, degree, reference_radius="STR(RADIUS)", is_internal=True)\n"\
"\n"\
"     By default when the 'is_internal' flag is set to True, evaluate\n"\
"     for the given 'radius' and 'reference_radius' relative radius power\n"\
"     series:\n"\
"       (reference_radius/radius)**(i+2) for i in range(0, degree+1) .\n"\
"\n"\
"     When the 'is_internal' flag is set to False, evaluate\n"\
"     for the given 'radius' and 'reference_radius' relative radius power\n"\
"     series:\n"\
"       (radius/reference_radius)**(i+1) for i in range(0, degree+1) .\n"\
"\n"


static PyObject* relradpow(PyObject *self, PyObject *args, PyObject *kwdict)
{
    static char *keywords[] = {
        "radius", "degree", "reference_radius", "is_internal", NULL
    };

    int degree;
    int is_internal;
    double rad, rad0 = RADIUS; // radius and reference radius
    PyObject *obj_is_internal = NULL; // boolean flag
    PyObject *arr_rrp = NULL; // P array
    PyObject *retval = NULL; // output tuple

    // parse input arguments
    if (!PyArg_ParseTupleAndKeywords(
        args, kwdict, "di|dO:relradpow", keywords,
        &rad, &degree, &rad0, &obj_is_internal
    ))
        goto exit;

    is_internal = (obj_is_internal == NULL) || PyObject_IsTrue(obj_is_internal);

    if (degree < 0)
    {
        PyErr_Format(PyExc_ValueError, "%s < 0", keywords[1]);
        goto exit;
    }

    if (rad < 0.0)
    {
        PyErr_Format(PyExc_ValueError, "%s < 0", keywords[0]);
        goto exit;
    }

    if (rad0 <= 0.0)
    {
        PyErr_Format(PyExc_ValueError, "%s <= 0", keywords[2]);
        goto exit;
    }

    // create the output array
    if (NULL == (arr_rrp = _get_new_double_array(1, NULL, degree+1)))
        goto exit;

    // evaluate the relative radius power series
    if (is_internal)
        shc_relradpow_internal(PyArray_DATA(arr_rrp), degree, rad/rad0);
    else
        shc_relradpow_external(PyArray_DATA(arr_rrp), degree, rad/rad0);

    retval = arr_rrp;

  exit:

    // decrease reference counters to the arrays
    if (!retval && arr_rrp){Py_DECREF(arr_rrp);}

    return retval;
}

#endif  /* PYMM_RELRADPOW_H */
