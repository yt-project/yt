"""
Tests for yt.funcs
"""
#-----------------------------------------------------------------------------
# Copyright (c) 2018, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------
from nose.tools import assert_raises

from yt import YTQuantity
from yt.funcs import validate_axis, validate_center
from yt.testing import assert_equal


def test_validate_axis():
    validate_axis(0)
    validate_axis('X')
    desired = ("Expected axis of int or char type ( can be "
               "[0, 1, 2, '0', '1', '2', 'x', 'y', 'z']), received ")

    with assert_raises(TypeError) as ex:
        validate_axis(3)
    assert_equal(str(ex.exception), desired + "'3'.")

    with assert_raises(TypeError) as ex:
        validate_axis('x-axis')
    assert_equal(str(ex.exception), desired + "'x-axis'.")

    with assert_raises(TypeError) as ex:
        validate_axis(['x', 'y'])
    assert_equal(str(ex.exception), desired + "'list'.")

def test_validate_center():
    validate_center("max")
    validate_center("MIN_")

    with assert_raises(TypeError) as ex:
        validate_center("avg")
    desired = ("Expected 'center' to be in ['c', 'center', 'm', 'max', 'min'] "
               "or the prefix to be 'max_'/'min_', received 'avg'.")
    assert_equal(str(ex.exception), desired)

    validate_center(YTQuantity(0.25, 'cm'))
    validate_center([0.25, 0.25, 0.25])

    class CustomCenter:
        def __init__(self, center):
            self.center = center

    with assert_raises(TypeError) as ex:
        validate_center(CustomCenter(10))
    desired = ("Expected 'center' to be a numeric object of type "
               "list/tuple/np.ndarray/YTArray/YTQuantity, received "
               "'yt.tests.test_funcs.test_validate_center.<locals>."
               "CustomCenter'.")
    assert_equal(str(ex.exception), desired)
