/*-----------------------------------------------------------------------------
 *
 * Magnetic Quasi Dipole Coordinates - C python bindings
 * -  sub-solar point evaluation
 *
 * Author: Martin Paces <martin.paces@eox.at>
 *
 *-----------------------------------------------------------------------------
 * Copyright (C) 2016 EOX IT Services GmbH
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

#ifndef PYQD_EVAL_SUBSOL_H
#define PYQD_EVAL_SUBSOL_H

#include "pymm_aux.h"
#include "qdipole/cqdipole.h"

/* python function definition */

#define DOC_EVAL_SUBSOL "\n"\
"   gdlat, gdlon = eval_subsol(time)\n"\
"     Inputs:\n"\
"       time  - MJD2000 time(s)\n"\
"     Outputs:\n"\
"       gdlat - sub-solar point latitude(s).\n"\
"       gdlon - sub-solar point longitudes(s).\n"\
""

static PyObject* eval_subsol(PyObject *self, PyObject *args, PyObject *kwdict)
{
    static char *keywords[] = {"time", NULL};

    PyObject *obj_time = NULL; // time object
    PyObject *arr_time = NULL; // time array
    PyObject *arr_gdlat = NULL; // gdlat array
    PyObject *arr_gdlon = NULL; // gdlon array
    PyObject *retval = NULL;

    // parse input arguments
    if (!PyArg_ParseTupleAndKeywords(
        args, kwdict, "O:eval_subsol", keywords, &obj_time
    ))
        goto exit;

    #define NPY_REQ (NPY_ALIGNED|NPY_CONTIGUOUS)

    // cast the objects to arrays
    if (NULL == (arr_time=_get_as_double_array(obj_time, 0, 1, NPY_REQ, keywords[0])))
        goto exit;

    // create the output arrays
    npy_intp ndim = PyArray_NDIM(arr_time);
    npy_intp *dims = PyArray_DIMS(arr_time);

    if (NULL == (arr_gdlat = PyArray_EMPTY(ndim, dims, NPY_DOUBLE, 0)))
        goto exit;

    if (NULL == (arr_gdlon = PyArray_EMPTY(ndim, dims, NPY_DOUBLE, 0)))
        goto exit;

    // evaluate the output values
    c_eval_subsol(
        (double*)PyArray_DATA(arr_gdlat),
        (double*)PyArray_DATA(arr_gdlon),
        (double*)PyArray_DATA(arr_time),
        ndim == 0 ? 1 : dims[0]
    );

    if (NULL == (retval = Py_BuildValue("NN", arr_gdlat, arr_gdlon))) goto exit;

  exit:

    // decrease reference counters to the arrays
    if (arr_time) {Py_DECREF(arr_time);}
    if (!retval)
    {
        if (arr_gdlat) {Py_DECREF(arr_gdlat);}
        if (arr_gdlon) {Py_DECREF(arr_gdlon);}
    }

    return retval;
}

#endif  /* PYQD_EVAL_SUBSOL_H */


