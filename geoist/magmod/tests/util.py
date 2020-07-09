#-------------------------------------------------------------------------------
#
#  Utility module tests.
#
# Author: Steve Shi Chen <chenshi80@gmail.com>
# 
# Original Author: Martin Paces <martin.paces@eox.at>
#-------------------------------------------------------------------------------
# Copyright (C) 2019 Geoist team
#
#-------------------------------------------------------------------------------
# pylint: disable=missing-docstring

from unittest import main, TestCase
from datetime import date, datetime
from numpy import array
from numpy.testing import assert_allclose
from geoist.magmod.util import vnorm, vincdecnorm, datetime_to_decimal_year


class FunctionTestMixIn(object):
    """ Simple function test mix in. """

    @staticmethod
    def eval(input_):
        raise NotImplementedError

    @staticmethod
    def _assert(tested_output, expected_output):
        assert_allclose(tested_output, expected_output, atol=1e-8)

    def test_accepted(self):
        for input_, expected_output in self.ACCEPTED:
            try:
                self._assert(self.eval(input_), expected_output)
            except AssertionError as exc:
                raise AssertionError("\ninput: %s\n%s" % (input_, exc))

    def test_rejected(self):
        for input_, expected_exception in self.REJECTED:
            try:
                self.assertRaises(expected_exception, self.eval, input_)
            except AssertionError as exc:
                raise AssertionError("\nInput: %s\n%s" % (input_, exc))


class TestDatetimeToDecimalYear(FunctionTestMixIn, TestCase):
    NAME = "datetime_to_decimal_year"

    @staticmethod
    def eval(input_):
        return datetime_to_decimal_year(input_)

    ACCEPTED = [
        (datetime(2001, 1, 12), 2001.0301369863014),
        (datetime(2012, 8, 31), 2012.6639344262296),
        (datetime(2014, 8, 31), 2014.66301369863),
        (datetime(2024, 12, 31, 23, 59, 59, 999), 2024.9999999684085),
    ]

    REJECTED = [
        (None, TypeError),
        (date(2001, 1, 12), TypeError),
    ]


class TestVnorm(FunctionTestMixIn, TestCase):
    NAME = "vnorm"

    @staticmethod
    def eval(input_):
        return vnorm(input_)

    ACCEPTED = [
        (1.0, 1.0),
        ([4, 3], 5),
        ([0, 0, 0], 0),
        (array([1, 0, 0]), 1),
        (array([0, 1.0, 0]), 1),
        (array([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]), array([1., 1., 1.])),
        (array([
            [
                [0.82544278, 0.23376316, 0.33205357],
                [0.06208131, 0.63504005, -0.78131731],
                [0.86222978, 0.65663201, 0.65170281],
                [-0.16437946, 0.40473597, -0.35760514],
            ],
            [
                [-0.48720026, 0.18952629, -0.8906141],
                [-0.06893185, -0.61986417, -0.79941557],
                [0.12094749, 0.33802827, 0.64999703],
                [-0.36555326, -0.1363087, 0.50091617],
            ],
            [
                [-0.25049911, 0.14976648, 0.80102138],
                [-0.96197621, 0.2603063, 0.17388555],
                [-0.94671758, -0.82102239, 0.62979664],
                [-0.99625197, -0.87799085, 0.95585576],
            ],
            [
                [0.04618681, -0.84110911, 0.97818344],
                [-0.01086467, 0.21928006, -0.96179247],
                [0.7898758, -0.86362723, 0.10181061],
                [0.14741587, -0.19785671, -0.06904178],
            ],
        ]), array([
            [0.91992422, 1.00875502, 1.26464317, 0.56454694],
            [1.03270411, 1.01392724, 0.74255474, 0.63492224],
            [0.85253449, 1.01162927, 1.40249626, 1.63616813],
            [1.29090689, 0.98653258, 1.17478559, 0.25621375],
        ]))
    ]
    REJECTED = [
        (None, TypeError),
    ]


class TestVincdecnorm(FunctionTestMixIn, TestCase):

    @staticmethod
    def eval(input_):
        return vincdecnorm(input_)

    @staticmethod
    def _assert(tested_output, expected_output):
        inclination, declination, intensity = tested_output
        inclination_ref, declination_ref, intensity_ref = expected_output
        try:
            assert_allclose(inclination, inclination_ref, atol=1e-8)
            assert_allclose(declination, declination_ref, atol=1e-8)
            assert_allclose(intensity, intensity_ref, atol=1e-8)
        except AssertionError as exc:
            template = "%s: received %s expected %s\n"
            raise AssertionError(
                (template % ("inclination", inclination, inclination_ref)) +
                (template % ("declination", declination, declination_ref)) +
                (template % ("intensity", intensity, intensity_ref)) +
                "%s" % exc
            )


    ACCEPTED = [
        ([1, 0, 0], (0, 0, 1)),
        ((0, 1, 0), (0, 90, 1)),
        (array([0, 0, 1]), (-90, 0, 1)),
        (
            array([[-1, 0, 0], [0, -1, 0], [0, 0, -1]]),
            (
                array([0, 0, 90]),
                array([180, -90, 0]),
                array([1, 1, 1]),
            )
        ),
        (
            array([
                [
                    [0.82544278, 0.23376316, 0.33205357],
                    [0.06208131, 0.63504005, -0.78131731],
                    [0.86222978, 0.65663201, 0.65170281],
                    [-0.16437946, 0.40473597, -0.35760514],
                ],
                [
                    [-0.48720026, 0.18952629, -0.8906141],
                    [-0.06893185, -0.61986417, -0.79941557],
                    [0.12094749, 0.33802827, 0.64999703],
                    [-0.36555326, -0.1363087, 0.50091617],
                ],
                [
                    [-0.25049911, 0.14976648, 0.80102138],
                    [-0.96197621, 0.2603063, 0.17388555],
                    [-0.94671758, -0.82102239, 0.62979664],
                    [-0.99625197, -0.87799085, 0.95585576],
                ],
                [
                    [0.04618681, -0.84110911, 0.97818344],
                    [-0.01086467, 0.21928006, -0.96179247],
                    [0.7898758, -0.86362723, 0.10181061],
                    [0.14741587, -0.19785671, -0.06904178],
                ],
            ]), (
                array([
                    [-21.15901257, 50.76300438, -31.01921101, 39.30418462],
                    [59.58823632, 52.03948549, -61.08670218, -52.0866531],
                    [-69.98055814, -9.89753007, -26.6830097, -35.74676721],
                    [-49.26615694, 77.14137388, -4.97166886, 15.63269911]
                ]),
                array([
                    [15.81197833, 84.41652502, 37.29112409, 112.10403758],
                    [158.74341437, -96.34549228, 70.3126441, -159.55035541],
                    [149.12596586, 164.85863574, -139.06717246, -138.61046712],
                    [-86.85694253, 92.83651461, -47.55387908, -53.31153904],
                ]),
                array([
                    [0.91992422, 1.00875502, 1.26464317, 0.56454694],
                    [1.03270411, 1.01392724, 0.74255474, 0.63492224],
                    [0.85253449, 1.01162927, 1.40249626, 1.63616813],
                    [1.29090689, 0.98653258, 1.17478559, 0.25621375],
                ])
            )
        )
    ]
    REJECTED = [
        (None, ValueError),
        (1.0, ValueError),
        ([1, 2], ValueError),
        ([1, 2, 3, 4], ValueError),
        (array([[1, 2]]), ValueError),
        (array([[1, 2, 3, 4]]), ValueError),
    ]

if __name__ == '__main__':
    main()
