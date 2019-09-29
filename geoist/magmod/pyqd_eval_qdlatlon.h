/*-----------------------------------------------------------------------------
 *
 * Magnetic Quasi Dipole Coordinates - C python bindings
 * -  QD coordinates' evaluation
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

#ifndef PYQD_EVAL_QDLATLON_H
#define PYQD_EVAL_QDLATLON_H

#include "pymm_aux.h"
#include "qdipole/cqdipole.h"

/* python function definition */

#define DOC_EVAL_QDLATLON "\n"\
"   qdlat, qdlon = eval_qdlatlon(gclat, gclon, gcrad, time, fname)\n"\
"   qdlat, qdlon, f11, f12, f21, f22, f__ = "\
                       "eval_qdlatlon(gclat, gclon, gcrad, time fname, True)\n"\
"     Inputs:\n"\
"       gclat - geocentric latitude(s).\n"\
"       gclon - geocentric longitudes(s).\n"\
"       gcrad - geocentric radial coordinate(s) in km.\n"\
"       time  - decimal year time(s)\n"\
"       fname - file-name of the model text file.\n"\
"       eval_base - boolean flag indicating whether the base vectors should\n"\
"                   be calculated or not (False by default).\n"\
"     Outputs:\n"\
"       qdlat - quasi-dipole latitude(s).\n"\
"       qdlon - quasi-dipole longitudes(s).\n"\
"       f11, f12, f21, f22 - base vectors.\n"\
"       f__    - F = |F1 x F2|\n"\
""

static PyObject* eval_qdlatlon(PyObject *self, PyObject *args, PyObject *kwdict)
{
    static char *keywords[] = {
        "gclat", "gclon", "gcrad", "time", "fname", "eval_base", NULL
    };

    int eval_base;
    const char *model_fname = NULL;
    PyObject *obj_gclat = NULL; // gclat object
    PyObject *obj_gclon = NULL; // gclon object
    PyObject *obj_gcrad = NULL; // gcrad object
    PyObject *obj_time = NULL; // time object
    PyObject *obj_flag = NULL; //  boolean flag object

    PyObject *arr_gclat = NULL; // gclat array
    PyObject *arr_gclon = NULL; // gclon array
    PyObject *arr_gcrad = NULL; // gcrad array
    PyObject *arr_time = NULL; // time array

    PyObject *arr_qdlat = NULL; // qdlat array
    PyObject *arr_qdlon = NULL; // qdlon array
    PyObject *arr_f11 = NULL; // f11 array
    PyObject *arr_f12 = NULL; // f12 array
    PyObject *arr_f21 = NULL; // f21 array
    PyObject *arr_f22 = NULL; // f22 array
    PyObject *arr_f = NULL; // f array
    PyObject *retval = NULL;

    // parse input arguments
    if (!PyArg_ParseTupleAndKeywords(
        args, kwdict, "OOOOs|O:eval_qdlatlon", keywords,
        &obj_gclat, &obj_gclon, &obj_gcrad, &obj_time, &model_fname, &obj_flag
    ))
        goto exit;

    eval_base = obj_flag ? (1 == PyObject_IsTrue(obj_flag)) : 0;

    #define NPY_REQ (NPY_ALIGNED|NPY_CONTIGUOUS)

    // cast the objects to arrays
    if (NULL == (arr_gclat=_get_as_double_array(obj_gclat, 0, 1, NPY_REQ, keywords[0])))
        goto exit;

    if (NULL == (arr_gclon=_get_as_double_array(obj_gclon, 0, 1, NPY_REQ, keywords[1])))
        goto exit;

    if (NULL == (arr_gcrad=_get_as_double_array(obj_gcrad, 0, 1, NPY_REQ, keywords[2])))
        goto exit;

    if (NULL == (arr_time=_get_as_double_array(obj_time, 0, 1, NPY_REQ, keywords[3])))
        goto exit;

    // check the dimensions
    npy_intp ndim = PyArray_NDIM(arr_gclat);
    npy_intp *dims = PyArray_DIMS(arr_gclat);

    if (_check_arr_dims_all_eq(arr_gclon, ndim, dims, keywords[1]))
        goto exit;

    if (_check_arr_dims_all_eq(arr_gcrad, ndim, dims, keywords[2]))
        goto exit;

    if (_check_arr_dims_all_eq(arr_time, ndim, dims, keywords[3]))
        goto exit;

    // create the output arrays
    if (NULL == (arr_qdlat = PyArray_EMPTY(ndim, dims, NPY_DOUBLE, 0)))
        goto exit;

    if (NULL == (arr_qdlon = PyArray_EMPTY(ndim, dims, NPY_DOUBLE, 0)))
        goto exit;

    if (eval_base)
    {
        if (NULL == (arr_f11 = PyArray_EMPTY(ndim, dims, NPY_DOUBLE, 0)))
            goto exit;

        if (NULL == (arr_f12 = PyArray_EMPTY(ndim, dims, NPY_DOUBLE, 0)))
            goto exit;

        if (NULL == (arr_f21 = PyArray_EMPTY(ndim, dims, NPY_DOUBLE, 0)))
            goto exit;

        if (NULL == (arr_f22 = PyArray_EMPTY(ndim, dims, NPY_DOUBLE, 0)))
            goto exit;

        if (NULL == (arr_f = PyArray_EMPTY(ndim, dims, NPY_DOUBLE, 0)))
            goto exit;
    }

    // evaluate the output values
    if (eval_base)
    {
        c_eval_qdlatlonvb(
            (double*)PyArray_DATA(arr_qdlat),
            (double*)PyArray_DATA(arr_qdlon),
            (double*)PyArray_DATA(arr_f11),
            (double*)PyArray_DATA(arr_f12),
            (double*)PyArray_DATA(arr_f21),
            (double*)PyArray_DATA(arr_f22),
            (double*)PyArray_DATA(arr_f),
            (double*)PyArray_DATA(arr_time),
            (double*)PyArray_DATA(arr_gcrad),
            (double*)PyArray_DATA(arr_gclat),
            (double*)PyArray_DATA(arr_gclon),
            ndim == 0 ? 1 : dims[0],
            model_fname
        );
    }
    else
    {
        c_eval_qdlatlon(
            (double*)PyArray_DATA(arr_qdlat),
            (double*)PyArray_DATA(arr_qdlon),
            (double*)PyArray_DATA(arr_time),
            (double*)PyArray_DATA(arr_gcrad),
            (double*)PyArray_DATA(arr_gclat),
            (double*)PyArray_DATA(arr_gclon),
            ndim == 0 ? 1 : dims[0],
            model_fname
        );
    }

    if (eval_base)
    {
        retval = Py_BuildValue(
            "NNNNNNN", arr_qdlat, arr_qdlon,
            arr_f11, arr_f12, arr_f21, arr_f22, arr_f
        );
    }
    else
    {
        retval = Py_BuildValue("NN", arr_qdlat, arr_qdlon);
    }

    if (NULL == retval) goto exit;

  exit:

    // decrease reference counters to the arrays
    if (arr_gclat) {Py_DECREF(arr_gclat);}
    if (arr_gclon) {Py_DECREF(arr_gclon);}
    if (arr_gcrad) {Py_DECREF(arr_gcrad);}
    if (arr_time) {Py_DECREF(arr_time);}
    if (!retval)
    {
        if (arr_qdlat) {Py_DECREF(arr_qdlat);}
        if (arr_qdlon) {Py_DECREF(arr_qdlon);}
        if (arr_f11) {Py_DECREF(arr_f11);}
        if (arr_f12) {Py_DECREF(arr_f12);}
        if (arr_f21) {Py_DECREF(arr_f21);}
        if (arr_f22) {Py_DECREF(arr_f22);}
        if (arr_f) {Py_DECREF(arr_f);}
    }

    return retval;
}

#endif  /* PYQD_EVAL_QDLATLON_H */

