#-------------------------------------------------------------------------------
#
#  Trace the field lines.
#
# Author: Martin Paces <martin.paces@eox.at>
#
#-------------------------------------------------------------------------------
# Copyright (C) 2018 EOX IT Services GmbH
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies of this Software or works derived from this Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#-------------------------------------------------------------------------------

try:
    # Python 2 - range as an iterator
    from buildins import xrange as range
except ImportError:
    pass

from numpy import asarray, copy
from .._pymm import GEOCENTRIC_CARTESIAN, GEOCENTRIC_SPHERICAL, convert
from ..util import vnorm, vrotate

EARTH_RADIUS = 6371.2
DEFAULT_MIN_RADIUS = 0.9977 * EARTH_RADIUS
DEFAULT_MAX_RADIUS = 1e6 * EARTH_RADIUS
DEFAULT_MAX_STEPS = 500
DEFAULT_STEP = 1e2


def trace_field_line(model, time, point,
                     input_coordinate_system=GEOCENTRIC_SPHERICAL,
                     output_coordinate_system=GEOCENTRIC_SPHERICAL,
                     trace_options=None, model_options=None):
    """ Trace field line passing trough the given point.
    Output:
        coords - array of field line Cartesian coordinates.
        field - array of field vectors in the requested in the Cartesian
                coordinate frame.

    Parameters:
        time -  MJD2000 time used by the field model.
        point - stating point (x, y, z) in kilometres.
        trace_options - options used by the field line tracing.
        model_options - dictionary of model options passed to the model.eval()
                method.
    """
    trace_options = trace_options or {}
    model_options = model_options or {}

    # convert the input point coordinate
    point = convert(point, input_coordinate_system, GEOCENTRIC_CARTESIAN)

    # trace the 'downstream' and 'upstream' segments
    coords_down, field_down = _trace_field_line(
        _model_function(model, time, **model_options),
        _convergence_condition(**trace_options), point, **trace_options
    )
    coords_up, field_up = _trace_field_line(
        _model_function(model, time, backward=True, **model_options),
        _convergence_condition(**trace_options), point, **trace_options
    )

    # join the north and reversed south segments
    coords_up.reverse()
    field_up.reverse()
    coords = asarray(coords_up[:-1] + coords_down)
    field = asarray(field_up[:-1] + field_down)

    # convert the output point coordinates and vector frame
    coords = convert(coords, GEOCENTRIC_CARTESIAN, output_coordinate_system)
    field = vrotate(
        field, None, coords, GEOCENTRIC_CARTESIAN, output_coordinate_system
    )

    return coords, field


def _model_function(model, time, backward=False, radius=EARTH_RADIUS,
                    **model_options):
    model_options["input_coordinate_system"] = GEOCENTRIC_CARTESIAN
    model_options["output_coordinate_system"] = GEOCENTRIC_CARTESIAN

    direction = -1.0 if backward else 1.0

    def _eval_model(point):
        field_vector = model.eval(time, point, **model_options)
        field_intensity = vnorm(field_vector)
        if field_intensity > 0:
            scale = direction * vnorm(point) / (radius *field_intensity)
        else:
            scale = 0
        return scale * field_vector, field_vector

    return _eval_model


def _convergence_condition(min_radius=DEFAULT_MIN_RADIUS,
                           max_radius=DEFAULT_MAX_RADIUS, **_):

    def _did_converge(point):
        radius = vnorm(point)
        return radius <= min_radius or radius >= max_radius

    return _did_converge


def _trace_field_line(eval_model, did_converge, start_point,
                      step=DEFAULT_STEP, max_steps=DEFAULT_MAX_STEPS, **_):
    """ Low level line tracing. """
    vectors = []
    points = [tuple(start_point)]
    point = copy(start_point)

    for _ in range(max_steps):
        direction, vector = eval_model(point)
        vectors.append(tuple(vector))

        point += step * direction
        points.append(tuple(point))

        if did_converge(point):
            break

    _, vector = eval_model(point)
    vectors.append(tuple(vector))

    return points, vectors
