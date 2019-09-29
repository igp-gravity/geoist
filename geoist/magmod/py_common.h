/*-----------------------------------------------------------------------------
 *
 *  common shared Python definitions Python 2/3 compatibility
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

#ifndef PY_COMMON
#define PY_COMMON

#include <Python.h>

/* checking the Python version */

#if PY_MAJOR_VERSION == 2
  #define IS_PYTHON_2
  #if PY_MINOR_VERSION < 6
    #error "Non-supported Python minor version!"
  #endif
#elif PY_MAJOR_VERSION == 3
  #define IS_PYTHON_3
  #if PY_MINOR_VERSION < 4
    #error "Non-supported Python minor version!"
  #endif
#else
  #error "Non-supported Python major version!"
#endif

/* Python 2/3 compatibility */

#ifdef IS_PYTHON_3
#define PyString_FromString PyUnicode_FromString
#endif

/* Python dictionary operations */

static void set_dict_item_str(PyObject *dict, const char * key, PyObject *value)
{
    PyDict_SetItemString(dict, key, value);
    Py_DECREF(value);
}

static void set_dict_item_str_long(PyObject *dict, const char * key, long value)
{
    set_dict_item_str(dict, key, PyLong_FromLong(value));
}

static void set_dict_item_str_double(PyObject *dict, const char * key, double value)
{
    set_dict_item_str(dict, key, PyFloat_FromDouble(value));
}

static void set_dict_item_str_str(PyObject *dict, const char * key, const char *value)
{
    #ifdef IS_PYTHON_2
    set_dict_item_str(dict, key, PyString_FromString(value));
    #else
    set_dict_item_str(dict, key, PyUnicode_FromString(value));
    #endif
}

/* Python module initialization */

#ifdef IS_PYTHON_3
static struct PyModuleDef module_definition = {
    PyModuleDef_HEAD_INIT, NULL, NULL, -1, NULL, NULL, NULL, NULL, NULL
};
#endif

PyObject* init_python_module(
    const char *name,
    const char *doc,
    PyMethodDef *methods
)
{
    PyObject *module = NULL;
#ifdef IS_PYTHON_3
    module_definition.m_name = name;
    module_definition.m_doc = doc;
    module_definition.m_methods = methods;
    module = PyModule_Create(&module_definition);
#else /* IS_PYTHON2 */
    module = Py_InitModule3(name, methods, doc);
#endif
    return module;
}

#endif /* PY_COMMON */
