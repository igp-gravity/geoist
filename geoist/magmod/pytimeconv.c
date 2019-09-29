/*-----------------------------------------------------------------------------
 *
 * Time conversion utilities.
 *
 * Author: Martin Paces <martin.paces@eox.at>
 *
 *-----------------------------------------------------------------------------
 * Copyright (C) 2018 EOX IT Services GmbH
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

#include "common.h" /* common definitions - to be included before Python.h */
#include <Python.h>
#include <numpy/arrayobject.h>

/* module version */
#include "version.h"

/* common python utilities */
#include "py_common.h"

/* time conversion subroutines */
#include "pytimeconv_mjd2000_to_decimal_year.h"
#include "pytimeconv_mjd2000_to_year_fraction.h"
#include "pytimeconv_decimal_year_to_mjd2000.h"

/*---------------------------------------------------------------------------*/
/* module's doc string */

#define DOC_PYTIMECONV \
"Time conversion utilities."

/*---------------------------------------------------------------------------*/
/*define module's methods */
static PyMethodDef pytimeconv_methods[] =
{
    {"mjd2000_to_decimal_year", (PyCFunction)pytimeconv_mjd2000_to_decimal_year, METH_VARARGS|METH_KEYWORDS, DOC_MJD2000_TO_DECIMAL_YEAR},
    {"mjd2000_to_year_fraction", (PyCFunction)pytimeconv_mjd2000_to_year_fraction, METH_VARARGS|METH_KEYWORDS, DOC_MJD2000_TO_YEAR_FRACTION},
    {"decimal_year_to_mjd2000", (PyCFunction)pytimeconv_decimal_year_to_mjd2000, METH_VARARGS|METH_KEYWORDS, DOC_DECIMAL_YEAR_TO_MJD2000},
    {NULL, NULL, 0, NULL} /* Sentinel - DO NOT REMOVE! */
} ;

/*---------------------------------------------------------------------------*/

/* module initialization  */
static PyObject* init_module(void)
{
    PyObject *module = init_python_module("_pytimeconv", DOC_PYTIMECONV, pytimeconv_methods);
    if (NULL == module)
        goto exit;

    PyObject *dict = PyModule_GetDict(module);
    if (NULL == dict)
        goto exit;

    /* metadata */
    set_dict_item_str_str(dict, "__author__", "Martin Paces (martin.paces@eox.at)");
    set_dict_item_str_str(dict, "__copyright__", "Copyright (C) 2018 EOX IT Services GmbH");
    set_dict_item_str_str(dict, "__licence__", "EOX licence (MIT style)");
    set_dict_item_str_str(dict, "__version__", VERSION);

  exit:
    return module;
}

#if PY_MAJOR_VERSION == 2

PyMODINIT_FUNC init_pytimeconv(void)
{
    import_array();
    init_module();
}

#else

PyObject* PyInit__pytimeconv(void)
{
    import_array();
    return init_module();
}

#endif
