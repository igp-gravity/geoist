/*-----------------------------------------------------------------------------
 *
 * Magnetic Quasi Dipole Coordinates - C python bindings
 *
 * Author: Martin Paces <martin.paces@eox.at>
 *
 *-----------------------------------------------------------------------------
 * Copyright (C) 2015 EOX IT Services GmbH
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies of this Software or works derived from this Software.
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

/* quasi dipole coordinates evaluation - old API with wrong MLT */
#include "pyqd_eval_apex.h"

/* Quasi-Dipole coordinates evaluation */
#include "pyqd_eval_qdlatlon.h"

/* Magnetic Local Time coordinates evaluation */
#include "pyqd_eval_mlt.h"

/* sub-solar point coordinates evaluation */
#include "pyqd_eval_subsol.h"

/*---------------------------------------------------------------------------*/
/* module's doc string */

#define DOC_PYQD \
"library evaluation magnetic quasi-dipole coordinates and local magnetic time."

/*---------------------------------------------------------------------------*/
/*define module's methods */
static PyMethodDef pyqd_methods[] =
{
    {"eval_qdlatlon", (PyCFunction)eval_qdlatlon, METH_VARARGS|METH_KEYWORDS, DOC_EVAL_QDLATLON},
    {"eval_mlt", (PyCFunction)eval_mlt, METH_VARARGS|METH_KEYWORDS, DOC_EVAL_MLT},
    {"eval_subsol", (PyCFunction)eval_subsol, METH_VARARGS|METH_KEYWORDS, DOC_EVAL_SUBSOL},
    {"eval_apex", (PyCFunction)eval_apex, METH_VARARGS|METH_KEYWORDS, DOC_EVAL_APEX},
    {NULL, NULL, 0, NULL} /* Sentinel - DO NOT REMOVE! */
} ;

/*---------------------------------------------------------------------------*/
/* module initialization  */

static PyObject* init_module(void)
{
    PyObject *module = init_python_module("_pyqd", DOC_PYQD, pyqd_methods);
    if (NULL == module)
        goto exit;

    PyObject *dict = PyModule_GetDict(module);
    if (NULL == dict)
        goto exit;

    /* add module specific exception */
    //PyDict_SetItemString(dict, "QDError", PyErr_NewException("_pyqd.QDError", NULL, NULL));

    /* metadata */
    set_dict_item_str_str(dict, "__author__", "Martin Paces (martin.paces@eox.at)");
    set_dict_item_str_str(dict, "__copyright__", "Copyright (C) 2015 EOX IT Services GmbH");
    set_dict_item_str_str(dict, "__licence__", "EOX licence (MIT style)");
    set_dict_item_str_str(dict, "__version__", VERSION);

  exit:
    return module;
}

/*---------------------------------------------------------------------------*/

#if PY_MAJOR_VERSION == 2

PyMODINIT_FUNC init_pyqd(void)
{
    import_array();
    init_module();
}

#else

PyObject* PyInit__pyqd(void)
{
    import_array();
    return init_module();
}

#endif
