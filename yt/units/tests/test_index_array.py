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
    result = index*index

    assert_equal(result.units, ((u.km**2).units, (u.g**2).units, (u.s**2).units))
    assert_equal(result.ndview, vals**2)
    assert(type(result.units) is tuple)

    index = IndexArray(vals, input_units=[u.km, u.g, u.s])

    result = index * u.km

    assert_equal(
        result.units, ((u.km**2).units, (u.g*u.km).units, (u.s*u.km).units))
    assert_equal(result.ndview, vals)
    assert(type(result.units) is tuple)

    result = u.km * index

    assert_equal(
        result.units, ((u.km**2).units, (u.g*u.km).units, (u.s*u.km).units))
    assert_equal(result.ndview, vals)
    assert(type(result.units) is tuple)

    result = index * 2

    assert_equal(result.units, index.units)
    assert_equal(result.ndview, 2*vals)
    assert(type(result.units) is tuple)

    result = 2 * index

    assert_equal(result.units, index.units)
    assert_equal(result.ndview, 2*vals)
    assert(type(result.units) is tuple)

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
    ret3 = ret1
    ret4 = YTQuantity(3, 'km')
    ret5 = YTQuantity(38, 's')
    ret6 = IndexArray([15, 16, 17], [u.km, u.g, u.s])

    sl1 = index[:, 0]
    sl2 = index[0, :]
    sl3 = index[..., 0]
    sl4 = index[1, 0]
    sl5 = index[12, 2]
    sl6 = index[5]

    compare_slicing(ret1, sl1, Unit, YTArray)
    compare_slicing(ret2, sl2, tuple, IndexArray)
    compare_slicing(ret3, sl3, Unit, YTArray)
    compare_slicing(ret4, sl4, Unit, YTQuantity)
    compare_slicing(ret5, sl5, Unit, YTQuantity)
    compare_slicing(ret6, sl6, tuple, IndexArray)

    index = IndexArray([0, 1, 2], input_units=[u.km, u.g, u.s])
    compare_slicing(index[0], YTQuantity(0, 'km'), Unit, YTQuantity)
    compare_slicing(index[2], YTQuantity(2, 's'), Unit, YTQuantity)
    compare_slicing(index[:], index, tuple, IndexArray)
    compare_slicing(index[...], index, tuple, IndexArray)

def test_boolean_unary_ops():
    vals = np.arange(30)
    vals.shape = (10, 3)
    index = IndexArray(vals, [u.km, u.km, u.km])

    assert_equal(np.isnan(index), np.isnan(vals))
    assert(all([i is np.bool_(False) for i in np.isnan(index).ravel()]))

def test_str_and_repr():
    vals = np.arange(6)
    vals.shape = (2, 3)
    index = IndexArray(vals, input_units=[u.km, u.g, u.s])

    assert_equal(str(index), '[[ 0.  1.  2.]\n [ 3.  4.  5.]] (km, g, s)')
    assert_equal(
        repr(index),
        'IndexArray([[ 0.,  1.,  2.],\n       [ 3.,  4.,  5.]]) (km, g, s)')
