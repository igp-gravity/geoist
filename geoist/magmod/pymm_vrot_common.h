/*-----------------------------------------------------------------------------
 *
 * Geomagnetic Model - C python bindings - vector rotation - common
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

#ifndef PYMM_VROT_COMMON_H
#define PYMM_VROT_COMMON_H

#include "pymm_aux.h"

/* array check*/

static int _vrot_arr_check(
    PyObject *arr_ref, PyObject *arr_checked,
    const char *label_ref, const char *label_checked
)
{
    //if ((PyArray_NDIM(arr_ref) > 1)||(PyArray_NDIM(arr_checked) > 0))
    if (PyArray_NDIM(arr_checked) > 0)
    {
        int d;

        if (PyArray_NDIM(arr_ref) != PyArray_NDIM(arr_checked)+1)
        {
            PyErr_Format(PyExc_ValueError, "Shape mismatch between '%s' and "
                "'%s'!", label_ref, label_checked);
            return 1;
        }

        for (d = 0; d < (PyArray_NDIM(arr_ref)-1); ++d)
        {
            if (PyArray_DIM(arr_ref, d) != PyArray_DIM(arr_checked, d))
            {
                PyErr_Format(PyExc_ValueError, "Shape mismatch between '%s' "
                    "and '%s'!", label_ref, label_checked);
                return 1;
            }
        }
    }

    return 0;
}

#endif  /* PYMM_VROT_COMMON_H */
