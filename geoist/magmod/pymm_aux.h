/*-----------------------------------------------------------------------------
 *
 * Geomagnetic Model - C python bindings - auxiliary functions
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

#ifndef PYMM_AUX_H
#define PYMM_AUX_H
#ifndef M_PI
#define M_PI 3.14159265358979323846  /* pi */
#endif

#ifndef M_1_PI
#define M_1_PI 0.31830988618379067154  /* 1/pi */
#endif
/*
 * Check the input python object and convert it to a double precision NumPy
 * array ensuring the native byte-order.
 */
static PyObject* _get_as_int_array(PyObject *data, int dmin, int dmax,
                int reqs, const char *label)
{
    PyArray_Descr *dtype = PyArray_DescrFromType(NPY_INT32);
    PyObject *arr = PyArray_FromAny(data, dtype, dmin, dmax, reqs, NULL);
    /*
    if (NULL == arr)
        PyErr_Format(PyExc_ValueError, "Failed to cast %s to an array!", label);
    */
    return arr;
}

/*
 * Check the input python object and convert it to a double precision NumPy
 * array ensuring the native byte-order.
 */
static PyObject* _get_as_double_array(PyObject *data, int dmin, int dmax,
                int reqs, const char *label)
{
    PyArray_Descr *dtype = PyArray_DescrFromType(NPY_FLOAT64);
    PyObject *arr = PyArray_FromAny(data, dtype, dmin, dmax, reqs, NULL);
    /*
    if (NULL == arr)
        PyErr_Format(PyExc_ValueError, "Failed to cast %s to an array!", label);
    */
    return arr;
}

/*
 * Get new allocated NumPy array. The first (N-1) dimensions as read from
 * the array of dimensions (allowing easily set the same shape as the input
 * matrix). The last Nth dimension is overridden by the 'dim_last' value.
 */
static PyObject* _get_new_double_array(npy_intp ndim, const npy_intp *dims, npy_intp dim_last)
{
    npy_intp i;
    npy_intp dims_new[MAX_OUT_ARRAY_NDIM];

    if (ndim > MAX_OUT_ARRAY_NDIM)
    {
        PyErr_Format(PyExc_ValueError, "The output array dimension "\
            " %d exceeds the allowed maximum value %d!", (int)ndim, MAX_OUT_ARRAY_NDIM);
        return NULL;
    }

    for (i = 0; i < (ndim-1); ++i)
        dims_new[i] = dims[i];
    if (ndim >= 1)
        dims_new[ndim-1] = dim_last;

    return PyArray_EMPTY(ndim, dims_new, NPY_DOUBLE, 0);
}

/*
 * Check that the dimension of the array has the required value.
 */
static int _check_array_dim_eq(PyObject *arr, int dim, size_t size, const char *label)
{
    if (dim < 0)
        dim += PyArray_NDIM(arr);
    int rv = PyArray_DIM(arr, dim) != size;
    if (rv)
        PyErr_Format(PyExc_ValueError, "The dimension #%d of '%s'"\
            " %ld is not equal the allowed value %ld!", dim, label,
            (size_t)PyArray_DIM(arr, dim), size);
    return rv;
}

/*
 * Check that the dimension of the array is greater then or equal to the required value.
 */
static int _check_array_dim_le(PyObject *arr, int dim, size_t size, const char *label)
{
    if (dim < 0)
        dim += PyArray_NDIM(arr);
    int rv = PyArray_DIM(arr, dim) < size;
    if (rv)
        PyErr_Format(PyExc_ValueError, "The dimension #%d of '%s'"\
            " %ld is lower than the minimum allowed value %ld!", dim, label,
            (size_t)PyArray_DIM(arr, dim), size);
    return rv;
}

/*
 * Check that the array dimensions match the required values
 */

static int _check_arr_dims_all_eq(PyObject *arr, npy_intp ndim, const npy_intp *dims, const char *label)
{
    npy_intp dim;

    if (PyArray_NDIM(arr) != ndim)
    {
        PyErr_Format(PyExc_ValueError, "The number of dimensions of '%s'"\
            " %ld does not match the required value %ld!", label,
            (size_t)(PyArray_NDIM(arr)), (size_t)ndim);
        return 1;
    }

    for (dim = 0; dim < ndim; ++dim)
    {
        if (PyArray_DIM(arr, dim) != dims[dim])
        {
            PyErr_Format(PyExc_ValueError, "The dimensions #%ld of '%s'"\
            " %ld does not match the required value %ld!", (size_t)dim, label,
            (size_t)PyArray_DIM(arr, dim), (size_t)dims[dim]);
            return 1;
        }
    }

    return 0;
}

/*
 * Check that two arrays have the same shape
 */

static int _check_equal_shape(PyObject *arr0, PyObject *arr1, const char *label0, const char *label1)
{
    npy_intp dim;

    if (PyArray_NDIM(arr0) != PyArray_NDIM(arr1))
    {
        PyErr_Format(
            PyExc_ValueError, "Array dimension mismatch between '%s' and '%s'!",
            label0, label1
        );
        return 1;
    }

    for (dim = 0; dim < PyArray_NDIM(arr0); ++dim)
    {
        if (PyArray_DIM(arr0, dim) != PyArray_DIM(arr1, dim))
        {
            PyErr_Format(
                PyExc_ValueError, "Array shape mismatch between '%s' and '%s'!",
                label0, label1
            );
            return 1;
        }
    }

    return 0;
}

/*
 *  Extract 1D Numpy array or broadcast scalar to a 1D C array.
 */
static int _extract_1d_double_array(double *out, size_t size, PyObject *arr_in, const char *label)
{
    npy_intp ndim = PyArray_NDIM(arr_in);
    void *data = PyArray_DATA(arr_in); //void *data = PyArray_DATA(arr_in);

    if (ndim == 0)
    {
        size_t i;
        for (i = 0 ; i < size; ++i)
            out[i] = *((double*)data);
    }
    else if ((ndim == 1) && (PyArray_DIM(arr_in, 0) == size))
    {
        npy_intp stride = PyArray_STRIDE(arr_in, 0);
        size_t i;
        for (i = 0 ; i < size; ++i)
            out[i] = *((double*)((char*)data + i*stride)); //*((double*)(data + i*stride));
    }
    else
    {
        PyErr_Format(PyExc_ValueError, "Invalid %s dimension!", label);
        return 1;
    }

    return 0;
}
/*
 * Extraction of the lower dimensional parts of the arrays.
 */

typedef struct {
    void *data; //void *data; 
    npy_intp ndim;
    const npy_intp *dim;
    const npy_intp *stride;
} ARRAY_DATA;

static ARRAY_DATA _array_to_arrd(PyObject *arr)
{
    ARRAY_DATA arrd = {
        PyArray_DATA(arr),
        PyArray_NDIM(arr),
        PyArray_DIMS(arr),
        PyArray_STRIDES(arr)
    };
    return arrd;
}

static ARRAY_DATA _get_arrd_item(const ARRAY_DATA *arrd, npy_intp idx)
{
    if (arrd->ndim > 0) // extract sub-dimension
    {
        ARRAY_DATA arrd_sub = {
            (char *)arrd->data + idx*arrd->stride[0], //arrd->data + idx*arrd->stride[0],
            arrd->ndim - 1,
            arrd->dim + 1,
            arrd->stride + 1
        };
        return arrd_sub;
    }
    else // treat as scalar
    {
        return *arrd;
    }
}

#endif  /* PYMM_AUX_H */
