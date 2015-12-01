"""
Test ndarray subclass that handles indexing along dimensions with units.




"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import yt.units as u
import numpy as np

from yt.testing import assert_equal
from yt.units.index_array import IndexArray
from yt.units.yt_array import YTArray, YTQuantity
from yt.units.unit_object import Unit

def test_multiplication():
    vals = np.random.random((100, 3))
    index = IndexArray(vals, input_units=[u.km, u.g, u.s])
    index *= index

    assert_equal(index.units[0], (u.km * u.km).units)
    assert_equal(index.units[1], (u.g * u.g).units)
    assert_equal(index.units[2], (u.s * u.s).units)
    assert(type(index.units) is tuple)

    index = IndexArray(vals, input_units=[u.km, u.g, u.s])

    index *= u.km

    assert_equal(index.units[0], (u.km * u.km).units)
    assert_equal(index.units[1], (u.km * u.g).units)
    assert_equal(index.units[2], (u.km * u.s).units)
    assert(type(index.units) is tuple)

    index2 = index * 2

    assert_equal(index.units[0], index2.units[0])
    assert_equal(index.units[1], index2.units[1])
    assert_equal(index.units[2], index2.units[2])
    assert(type(index.units) is tuple)

    index3 = 2 * index

    assert_equal(index.units[0], index3.units[0])
    assert_equal(index.units[1], index3.units[1])
    assert_equal(index.units[2], index3.units[2])
    assert(type(index.units) is tuple)

def compare_slicing(desired, actual, unit_type, array_type):
    assert_equal(desired, actual)
    assert(type(actual.units) is unit_type)
    assert(type(actual) is array_type)

def test_slicing():
    vals = np.arange(300)
    vals.shape = (100, 3)
    index = IndexArray(vals, input_units=[u.km, u.g, u.s])

    ret1 = YTArray(3*np.arange(100), 'km')
    ret2 = IndexArray(np.arange(3), [u.km, u.g, u.s])
    ret3 = YTQuantity(3, 'km')
    ret4 = YTQuantity(38, 's')

    sl1 = index[:, 0]
    sl2 = index[0, :]
    sl3 = index[1, 0]
    sl4 = index[12, 2]

    compare_slicing(ret1, sl1, Unit, YTArray)
    compare_slicing(ret2, sl2, tuple, IndexArray)
    compare_slicing(ret3, sl3, Unit, YTQuantity)
    compare_slicing(ret4, sl4, Unit, YTQuantity)

def test_str_and_repr():
    vals = np.arange(6)
    vals.shape = (2, 3)
    index = IndexArray(vals, input_units=[u.km, u.g, u.s])

    assert_equal(str(index), '[[0 1 2]\n [3 4 5]] (km, g, s)')
    assert_equal(
        repr(index), 'IndexArray([[0, 1, 2],\n       [3, 4, 5]]) (km, g, s)')
