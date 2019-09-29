/** 
 * @file math_aux.h
 * @author Martin Paces <martin.paces@eox.at>
 * @brief Auxiliary mathematical subroutines.
 *
 * This file contains definitions of auxiliary mathematical subroutines 
 * to be used all over the whole program.
 *
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
 */

#ifndef MATH_AUX_H
#define MATH_AUX_H 1

#include<math.h>

/**
 * @brief Evaluate quadratic norm of a 3D vector.
 *  v = (x, y, z)
 *  l = |v| = sqrt(x^2 + y^2 + z^2) 
 */
static double norm3d(double x, double y, double z)
{ 
    return sqrt(x*x + y*y + z*z); 
}

/**
 * @brief Evaluate quadratic norm of a 2D vector.
 *  v = (x, y)
 *  l = |v| = sqrt(x^2 + y^2) 
 */
static double norm2d(double x, double y)
{
    return sqrt(x*x + y*y);
}

/**
 * @brief Evaluate derivative of a quadratic norm of a 3D vector.
 *  l = |v| 
 *  dl = d(|v|)
 *  v = (x, y, z)
 *  dv = (dx, dy, dz) 
 *  if |v| >= 0:
 *    dl = (v/|v|)*dv = (x*dx + y*dy + z*dz)/l
 *  if |v| == 0: 
 *    dl = |dv| = sqrt(dx^2 + dy^2 + dz^2)
 */
static double dnorm3d(double x, double y, double z, double l,
               double dx, double dy, double dz)
{ 
    if (l > 0.0)
        return (x*dx + y*dy + z*dz)/l;
    else if (l == 0.0)
        return norm3d(dx, dy, dz);
    else 
        return NAN;
}

/**
 * @brief Evaluate derivative of a quadratic norm of a 2D vector.
 *  l = |v| 
 *  dl = d(|v|)
 *  v = (x, y)
 *  dv = (dx, dy) 
 *  if |v| >= 0:
 *    dl = (v/|v|)*dv = (x*dx + y*dy)/l
 *  if |v| == 0: 
 *    dl = |dv| = sqrt(dx^2 + dy^2)
 */
static double dnorm2d(double x, double y, double l, double dx, double dy)
{ 
    if (l > 0.0)
        return (x*dx + y*dy)/l;
    else if (l == 0.0)
        return norm2d(dx, dy);
    else 
        return NAN;
}

/**
 * @brief 2D vector rotation.
 *
 * x' = x*cos(a) - y*sin(a)
 * y' = x*sin(a) + y*cos(a)
 */
static void rot2d(double *x_out, double *y_out, double x, double y, double a_sin, double a_cos)
{
    *x_out = x*a_cos - y*a_sin;
    *y_out = x*a_sin + y*a_cos;
}

#endif  /*MATH_AUX_H*/
