/*-----------------------------------------------------------------------------
 *
 * Geomagnetic Model - C python bindings
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

#include "common.h" /* common definitions - to be included before Python.h */
#include <Python.h>
#include <numpy/arrayobject.h>

/* module version */
#include "version.h"

/* common python utilities */
#include "py_common.h"

/* Coordinate conversions. */
#include "pymm_cconv.h"

/* spherical harmonic model evaluation */
#include "pymm_sheval.h"

/* evaluation of the associative Legendre functions */
#include "pymm_legendre.h"

/* evaluation of the relative radius powers */
#include "pymm_relradpow.h"

/* evaluation of the series of longitude sine and cosine values */
#include "pymm_lonsincos.h"

/* final spherical-harmonic gradient evaluation */
#include "pymm_sphargrd.h"

/* final spherical-harmonic gradient potential */
#include "pymm_spharpot.h"

/* vector rotations */
#include "pymm_vrot_sph2geod.h"
#include "pymm_vrot_sph2cart.h"
#include "pymm_vrot_cart2sph.h"

/*---------------------------------------------------------------------------*/
/* module's doc string */

#define DOC_PYMM \
"This module provides bindings to the Geomagnetic Model library."

/*---------------------------------------------------------------------------*/
/*define module's methods */
static PyMethodDef pymm_methods[] =
{
    {"vrot_sph2cart", (PyCFunction)vrot_sph2cart, METH_VARARGS|METH_KEYWORDS, DOC_VROT_SPH2CART},
    {"vrot_cart2sph", (PyCFunction)vrot_cart2sph, METH_VARARGS|METH_KEYWORDS, DOC_VROT_CART2SPH},
    {"vrot_sph2geod", (PyCFunction)vrot_sph2geod, METH_VARARGS|METH_KEYWORDS, DOC_VROT_SPH2GEOD},
    {"spharpot", (PyCFunction)spharpot, METH_VARARGS|METH_KEYWORDS, DOC_SPHARPOT},
    {"sphargrd", (PyCFunction)sphargrd, METH_VARARGS|METH_KEYWORDS, DOC_SPHARGRD},
    {"lonsincos", (PyCFunction)lonsincos, METH_VARARGS|METH_KEYWORDS, DOC_LONSINCOS},
    {"relradpow", (PyCFunction)relradpow, METH_VARARGS|METH_KEYWORDS, DOC_RELRADPOW},
    {"legendre", (PyCFunction)legendre, METH_VARARGS|METH_KEYWORDS, DOC_LEGENDRE},
    {"sheval", (PyCFunction)sheval, METH_VARARGS|METH_KEYWORDS, DOC_SHEVAL},
    {"convert", (PyCFunction)convert, METH_VARARGS|METH_KEYWORDS, DOC_CONVERT},
    {NULL, NULL, 0, NULL} /* Sentinel - DO NOT REMOVE! */
} ;

/*---------------------------------------------------------------------------*/
/* module initialization  */

static PyObject* init_module(void)
{
    PyObject *module = init_python_module("_pymm", DOC_PYMM, pymm_methods);
    if (NULL == module)
        goto exit;

    PyObject *dict = PyModule_GetDict(module);
    if (NULL == dict)
        goto exit;

    /* add module specific exception */
    //set_dict_item_str(dict, "MMError", PyErr_NewException("_pymm.MMError", NULL, NULL));

    /* integer constants */
    set_dict_item_str_long(dict, "GEODETIC_ABOVE_WGS84", CT_GEODETIC_ABOVE_WGS84);
    set_dict_item_str_long(dict, "GEOCENTRIC_SPHERICAL", CT_GEOCENTRIC_SPHERICAL);
    set_dict_item_str_long(dict, "GEOCENTRIC_CARTESIAN", CT_GEOCENTRIC_CARTESIAN);
    set_dict_item_str_long(dict, "POTENTIAL", SM_POTENTIAL);
    set_dict_item_str_long(dict, "GRADIENT", SM_GRADIENT);
    set_dict_item_str_long(dict, "POTENTIAL_AND_GRADIENT", SM_POTENTIAL_AND_GRADIENT);
    set_dict_item_str_double(dict, "EARTH_RADIUS", RADIUS);

    /* module metadata */
    set_dict_item_str_str(dict, "__author__", "Martin Paces (martin.paces@eox.at)");
    set_dict_item_str_str(dict, "__copyright__", "Copyright (C) 2014 EOX IT Services GmbH");
    set_dict_item_str_str(dict, "__licence__", "EOX licence (MIT style)");
    set_dict_item_str_str(dict, "__version__", VERSION);

  exit:
    return module;
}

#if PY_MAJOR_VERSION == 2

PyMODINIT_FUNC init_pymm(void)
{
    import_array();
    init_module();
}

#else

PyObject* PyInit__pymm(void)
{
    import_array();
    return init_module();
}

#endif
