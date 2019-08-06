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
from yt.testing import assert_equal, fake_amr_ds


def test_validate_axis():
    validate_axis(None, 0)
    validate_axis(None, 'X')

    ds = fake_amr_ds(geometry="cylindrical")
    ds.slice("Theta", 0.25)

    with assert_raises(TypeError) as ex:
        # default geometry is cartesian
        ds = fake_amr_ds()
        ds.slice("r", 0.25)
    desired = ("Expected axis of int or char type (can be "
               "[0, 'x', 'X', 1, 'y', 'Y', 2, 'z', 'Z']), received 'r'.")
    assert_equal(str(ex.exception)[:40], desired[:40])

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
    assert_equal(str(ex.exception)[:50], desired[:50])
