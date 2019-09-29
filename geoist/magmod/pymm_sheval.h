/*-----------------------------------------------------------------------------
 *
 * Geomagnetic Model - C python bindings - spherical harmonic model evaluation
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

#ifndef PYMM_SHEVAL_H
#define PYMM_SHEVAL_H

#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "math_aux.h"
#include "geo_conv.h"
#include "shc.h"
#include "pymm_aux.h"
#include "pymm_coord.h"
#include "pymm_cconv.h"

#ifndef NAN
#define NAN (0.0/0.0)
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846  /* pi */
#endif

#ifndef M_1_PI
#define M_1_PI 0.31830988618379067154  /* 1/pi */
#endif
/* Earth radius in km */
#ifndef RADIUS
#define RADIUS  6371.2
#endif

#ifndef STR
#define SRINGIFY(x) #x
#define STR(x) SRINGIFY(x)
#endif

/* mode of evaluation enumerate */
typedef enum {
    SM_INVALID = 0x0,
    SM_POTENTIAL = 0x1,
    SM_GRADIENT = 0x2,
    SM_POTENTIAL_AND_GRADIENT = 0x3,
} SHEVAL_MODE;

/* mode of evaluation check */
static SHEVAL_MODE _check_sheval_mode(int mode, const char *label)
{
    switch (mode)
    {
        case SM_POTENTIAL:
            return SM_POTENTIAL;
        case SM_GRADIENT:
            return SM_GRADIENT;
        case SM_POTENTIAL_AND_GRADIENT:
            return SM_POTENTIAL_AND_GRADIENT;
        default:
            PyErr_Format(PyExc_ValueError, "Invalid mode '%s'!", label);
            return SM_INVALID;
    }
}

/* magnetic model - auxiliary structure */
typedef struct {
    int is_internal;
    int degree;
    int nterm;
    int coord_in;
    int coord_out;
    double elps_a;
    double elps_eps2;
    double crad_ref;
    double clat_last;
    double clon_last;
    double crad_last;
    double scale_potential;
    double scale_gradient0;
    double scale_gradient1;
    double scale_gradient2;
    const double *cg;
    const double *ch;
    double *lp;
    double *ldp;
    double *lsin;
    double *lcos;
    double *rrp;
    double *psqrt;
} MODEL;

/* model structure destruction */
static void _sheval_model_destroy(MODEL *model)
{
    if(NULL != model->lp)
        free(model->lp);
    if(NULL != model->ldp)
        free(model->ldp);
    if(NULL != model->lsin)
        free(model->lsin);
    if(NULL != model->lcos)
        free(model->lcos);
    if(NULL != model->rrp)
        free(model->rrp);
    if(NULL != model->psqrt)
        free(model->psqrt);

    memset(model, 0, sizeof(MODEL));
}

/* model structure initialization */
static int _sheval_model_init(MODEL *model, int is_internal, int degree,
    int coord_in, int coord_out, const double *cg, const double *ch,
    const double scale_potential, const double *scale_gradient
    )
{
    const int nterm = ((degree+1)*(degree+2))/2;

    memset(model, 0, sizeof(MODEL));

    model->is_internal = is_internal;
    model->degree = degree;
    model->nterm = nterm;
    model->coord_in = coord_in;
    model->coord_out = coord_out;
    model->elps_a = WGS84_A;
    model->elps_eps2 = WGS84_EPS2;
    model->crad_ref = RADIUS;
    model->clat_last = NAN;
    model->clon_last = NAN;
    model->crad_last = NAN;
    model->scale_potential = scale_potential;
    model->scale_gradient0 = scale_gradient[0];
    model->scale_gradient1 = scale_gradient[1];
    model->scale_gradient2 = scale_gradient[2];
    model->cg = cg;
    model->ch = ch;

    if (NULL == (model->lp = (double*)malloc(nterm*sizeof(double))))
        goto error;

    if (NULL == (model->ldp = (double*)malloc(nterm*sizeof(double))))
        goto error;

    if (NULL == (model->lsin = (double*)malloc((degree+1)*sizeof(double))))
        goto error;

    if (NULL == (model->lcos = (double*)malloc((degree+1)*sizeof(double))))
        goto error;

    if (NULL == (model->rrp = (double*)malloc((degree+1)*sizeof(double))))
        goto error;

    if (NULL == (model->psqrt = shc_presqrt(degree)))
        goto error;

    return 0;

  error:
    _sheval_model_destroy(model);
    return 1;
}

/* single point evaluation */
static void _sheval_model_eval(MODEL *model, int mode, double *fpot,
    double *fx, double *fy, double *fz, double x, double y, double z)
{
    double glat, glon, ghgt;
    double clat, clon, crad;
    double flat, flon, frad;
    double tmp;

    int glat_needed = (model->coord_out == CT_GEODETIC_ABOVE_WGS84);

    void (*shc_relradpow)(double *rrp, int degree, double relrad) = (
        model->is_internal ? shc_relradpow_internal : shc_relradpow_external
    );

    // convert the input coordinates
    switch(model->coord_in)
    {
        case CT_GEODETIC_ABOVE_WGS84:
            glat = x; glon = y; ghgt = z;
            geodetic2geocentric_sph(&crad, &clat, &clon, glat, glon, ghgt, model->elps_a, model->elps_eps2);
            break;
        case CT_GEOCENTRIC_SPHERICAL:
            clat = DG2RAD*x; clon = DG2RAD*y; crad = z;
            if (glat_needed) {
                geocentric_sph2geodetic(&glat, &glon, &ghgt, crad, clat, clon, model->elps_a, model->elps_eps2);
            }
            break;
        case CT_GEOCENTRIC_CARTESIAN:
            cart2sph(&crad, &clat, &clon, x, y, z);
            if (glat_needed) {
                geocentric_sph2geodetic(&glat, &glon, &ghgt, crad, clat, clon, model->elps_a, model->elps_eps2);
            }
            break;
        default:
            return;
    }

    // associative Legendre functions
    if (model->clat_last != clat)
        shc_legendre(model->lp, model->ldp, model->degree, clat, model->psqrt);

    // relative radial powers
    if (model->clon_last != clon)
        shc_azmsincos(model->lsin, model->lcos, model->degree, clon);

    // sines/cosines of the longitude
    if (model->crad_last != crad)
        shc_relradpow(model->rrp, model->degree, crad/model->crad_ref);

    // save the last evaluated coordinate
    model->clat_last = clat;
    model->clon_last = clon;
    model->crad_last = crad;

    // evaluate the model
    shc_eval(
        fpot, &flat, &flon, &frad, model->degree, mode,
        clat, crad, model->cg, model->ch, model->lp, model->ldp,
        model->rrp, model->lsin, model->lcos, model->is_internal
    );

    // project the produced vectors to the desired coordinate system
    if (mode&0x2)
    {
        switch(model->coord_out)
        {
            case CT_GEOCENTRIC_SPHERICAL:
                *fx = flat;
                *fy = flon;
                *fz = frad;
                break;

            case CT_GEODETIC_ABOVE_WGS84:
                tmp = clat - DG2RAD*glat;
                rot2d(fz, fx, frad, flat, sin(tmp), cos(tmp));
                *fy = flon;
                break;

            case CT_GEOCENTRIC_CARTESIAN:
                rot2d(&tmp, fz, frad, flat, sin(clat), cos(clat));
                rot2d(fx, fy, tmp, flon, sin(clon), cos(clon));
                break;
        }

        // final scaling
        *fx *= model->scale_gradient0;
        *fy *= model->scale_gradient1;
        *fz *= model->scale_gradient2;
    }

    if (mode&0x1)
    {
        // final scaling
        *fpot *= model->scale_potential;
    }
}

/* recursive batch model_evaluation for the input numpy arrays of coordinates */
#define S(a) ((double*)(a).data)
#define P(a,i) ((double*)((char*)(a).data+(i)*(a).stride[0])) //((double*)((a).data+(i)*(a).stride[0]))
#define V(a,i) (*P(a,i))

static void _sheval1(ARRAY_DATA arrd_in, ARRAY_DATA arrd_pot, MODEL *model)
{
    if (arrd_in.ndim > 1)
    {
        npy_intp i;
        for(i = 0; i < arrd_in.dim[0]; ++i)
            _sheval1(_get_arrd_item(&arrd_in, i), _get_arrd_item(&arrd_pot, i), model);
        return;
    }

    _sheval_model_eval(model, SM_POTENTIAL, S(arrd_pot), NULL, NULL, NULL,
                    V(arrd_in, 0), V(arrd_in, 1), V(arrd_in, 2));
}

static void _sheval2(ARRAY_DATA arrd_in, ARRAY_DATA arrd_grd, MODEL *model)
{
    if (arrd_in.ndim > 1)
    {
        npy_intp i;
        for(i = 0; i < arrd_in.dim[0]; ++i)
            _sheval2(_get_arrd_item(&arrd_in, i), _get_arrd_item(&arrd_grd, i), model);
        return;
    }

    _sheval_model_eval(model, SM_GRADIENT, NULL, P(arrd_grd, 0), P(arrd_grd, 1),
        P(arrd_grd, 2), V(arrd_in, 0), V(arrd_in, 1), V(arrd_in, 2));

}

static void _sheval3(ARRAY_DATA arrd_in, ARRAY_DATA arrd_pot, ARRAY_DATA arrd_grd, MODEL *model)
{
    if (arrd_in.ndim > 1)
    {
        npy_intp i;
        for(i = 0; i < arrd_in.dim[0]; ++i)
            _sheval3(_get_arrd_item(&arrd_in, i), _get_arrd_item(&arrd_pot, i),
                    _get_arrd_item(&arrd_grd, i), model);
        return;
    }

    _sheval_model_eval(model, SM_POTENTIAL_AND_GRADIENT, S(arrd_pot),
        P(arrd_grd, 0), P(arrd_grd, 1), P(arrd_grd, 2), V(arrd_in, 0),
        V(arrd_in, 1), V(arrd_in, 2));
}

#undef V
#undef P
#undef S

/* python function definition */

#define DOC_SHEVAL "\n"\
"   arr_out = sheval(arr_in, degree, coef_g, coef_h, coord_type_in=GEODETIC_ABOVE_WGS84,\n"\
"                    coord_type_out=GEODETIC_ABOVE_WGS84, mode=GRADIENT, is_internal=True,\n"\
"                    scale_potential=1.0, scale_gradient=1.0)\n"\
"\n"\
"     Parameters:\n"\
"       arr_in - array of 3D coordinates (up to 16 dimensions).\n"\
"       degree - degree of the spherical harmonic model.\n"\
"       coef_g - vector of spherical harmonic model coefficients.\n"\
"       coef_h - vector of spherical harmonic model coefficients.\n"\
"       coord_type_in - type of the input coordinates.\n"\
"       mode - quantity to be evaluated:\n"\
"                  POTENTIAL\n"\
"                  GRADIENT (default)\n"\
"                  POTENTIAL_AND_GRADIENT\n"\
"       is_internal - boolean flag set to True by default. When set to False \n"\
"                     external field evaluation is used.\n"\
"       scale_potential - scalar value multiplied with the result potentials.\n"\
"       scale_gradient - scalar or 3 element array multiplied with the result\n"\
"                        gradient components.\n"\
"\n"


static PyObject* sheval(PyObject *self, PyObject *args, PyObject *kwdict)
{
    static char *keywords[] = {
        "arr_in", "degree", "coef_g", "coef_h", "coord_type_in",
        "coord_type_out", "mode", "is_internal",
        "scale_potential", "scale_gradient", NULL
    };
    int ct_in = CT_GEODETIC_ABOVE_WGS84;
    int ct_out = CT_GEODETIC_ABOVE_WGS84;
    int nterm, degree = 0, mode = 0x2, is_internal;
    double scale_potential = 1.0;
    double scale_gradient[3] = {1.0, 1.0, 1.0};
    PyObject *obj_is_internal = NULL; // boolean flag
    PyObject *obj_in = NULL; // input object
    PyObject *obj_cg = NULL; // coef_g object
    PyObject *obj_ch = NULL; // coef_h object
    PyObject *obj_scale = NULL; // gradient scale object
    PyObject *arr_in = NULL; // input array
    PyObject *arr_cg = NULL; // coef_g array
    PyObject *arr_ch = NULL; // coef_h array
    PyObject *arr_pot = NULL; // output array
    PyObject *arr_grd = NULL; // output array
    PyObject *arr_scale = NULL; // gradient scale array
    PyObject *retval = NULL;
    MODEL model;

    // parse input arguments
    if (!PyArg_ParseTupleAndKeywords(
        args, kwdict, "OiOO|iiiOdO:sheval", keywords,
        &obj_in, &degree, &obj_cg, &obj_ch, &ct_in, &ct_out, &mode,
        &obj_is_internal, &scale_potential, &obj_scale
    ))
        goto exit;

    is_internal = (obj_is_internal == NULL) || PyObject_IsTrue(obj_is_internal);

    // check the type of the coordinate transformation
    if (CT_INVALID == _check_coord_type(ct_in, keywords[4]))
        goto exit;

    if (CT_INVALID == _check_coord_type(ct_out, keywords[5]))
        goto exit;

    // check the operation mode
    if (SM_INVALID == _check_sheval_mode(mode, keywords[6]))
        goto exit;

    // cast the objects to arrays
    if (NULL == (arr_in=_get_as_double_array(obj_in, 1, 0, NPY_ALIGNED, keywords[0])))
        goto exit;

    if (NULL == (arr_cg=_get_as_double_array(obj_cg, 1, 1, NPY_IN_ARRAY, keywords[2])))
        goto exit;

    if (NULL == (arr_ch=_get_as_double_array(obj_ch, 1, 1, NPY_IN_ARRAY, keywords[3])))
        goto exit;

    if (degree < 0)
    {
        PyErr_Format(PyExc_ValueError, "Invalid value %d of '%s'!", degree, keywords[1]);
        goto exit;
    }

    // check dimensions of the coefficient arrays
    nterm = ((degree+1)*(degree+2))/2;

    // check maximum allowed input array dimension
    if (PyArray_NDIM(arr_in) > MAX_OUT_ARRAY_NDIM)
    {
        PyErr_Format(PyExc_ValueError, "The input array dimension of '%s'"\
            " %d exceeds the allowed maximum value %d!", keywords[0],
            PyArray_NDIM(arr_in), MAX_OUT_ARRAY_NDIM);
        goto exit;
    }

    // check the last dimension (required length of the coordinates vectors)
    if (_check_array_dim_eq(arr_in, -1, 3, keywords[0]))
        goto exit;

    // check the coefficients arrays' dimensions
    if (_check_array_dim_le(arr_cg, 0, nterm, keywords[2]))
        goto exit;

    if (_check_array_dim_le(arr_ch, 0, nterm, keywords[3]))
        goto exit;

    // handle the gradient scale factors
    if (NULL != obj_scale)
    {
        if (NULL == (arr_scale=_get_as_double_array(obj_scale, 0, 1, NPY_IN_ARRAY, keywords[9])))
            goto exit;

        if (_extract_1d_double_array(scale_gradient, 3, arr_scale, keywords[9]))
            goto exit;
    }

    // create the output arrays
    if (mode & SM_POTENTIAL)
        if (NULL == (arr_pot = PyArray_EMPTY(PyArray_NDIM(arr_in)-1, PyArray_DIMS(arr_in), NPY_DOUBLE, 0)))
            goto exit;

    if (mode & SM_GRADIENT)
        if (NULL == (arr_grd = _get_new_double_array(PyArray_NDIM(arr_in), PyArray_DIMS(arr_in), 3)))
            goto exit;

    // evaluate the model

    if(_sheval_model_init(
        &model, is_internal, degree, ct_in, ct_out,
        PyArray_DATA(arr_cg), PyArray_DATA(arr_ch),
        scale_potential, scale_gradient
    ))
        goto exit;

    switch(mode)
    {
        case SM_POTENTIAL:
            _sheval1(_array_to_arrd(arr_in), _array_to_arrd(arr_pot), &model);
            retval = arr_pot;
            break;

        case SM_GRADIENT:
            _sheval2(_array_to_arrd(arr_in), _array_to_arrd(arr_grd), &model);
            retval = arr_grd;
            break;

        case SM_POTENTIAL_AND_GRADIENT:
            _sheval3(_array_to_arrd(arr_in), _array_to_arrd(arr_pot), _array_to_arrd(arr_grd), &model);
            if (NULL == (retval = Py_BuildValue("NN", arr_pot, arr_grd)))
                goto exit;
            break;
    }

    _sheval_model_destroy(&model);

  exit:

    // decrease reference counters to the arrays
    if (arr_in){Py_DECREF(arr_in);}
    if (arr_cg){Py_DECREF(arr_cg);}
    if (arr_ch){Py_DECREF(arr_ch);}
    if (!retval && arr_grd){Py_DECREF(arr_grd);}
    if (!retval && arr_pot){Py_DECREF(arr_pot);}

    return retval;
}

#endif  /* PYMM_SHEVAL_H */
