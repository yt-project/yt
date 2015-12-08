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
import operator

from yt.testing import \
    assert_equal, \
    assert_almost_equal, \
    assert_raises, \
    fake_random_ds
from yt.units.index_array import IndexArray
from yt.units.yt_array import YTArray, YTQuantity
from yt.units.unit_object import Unit

def test_creation_from_ytarray():
    arr = YTArray([1, 2, 3], 'g')
    iarr = IndexArray(arr)

    assert_equal(iarr.units, (u.g.units, )*3)
    assert_equal(iarr.ndview, np.array([1, 2, 3]))
    assert(type(iarr.units) is tuple)

def test_creation_errors():
    assert_raises(RuntimeError, IndexArray, [1, 2, 3, 4], ['g', 'cm', 's'])
    assert_raises(NotImplementedError, IndexArray, np.random.random((10, 10, 3)),
                  'g')

def test_registry_association():
    ds = fake_random_ds(64)

    arr = IndexArray([1, 2, 3], ('g', 'km', 's'), registry=ds.unit_registry)

    assert(all([u.registry is ds.unit_registry for u in (3*arr).units]))
    assert(all([u.registry is ds.unit_registry for u in (arr*3).units]))

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

def test_addition():
    vals = np.random.random((100, 3))
    index = IndexArray(vals, input_units=[u.km, u.g, u.s])

    row_vals = np.random.random(3)
    row_index = IndexArray(row_vals, input_units=[u.km, u.g, u.s])

    ret1 = index + index

    assert_equal(ret1.units, index.units)
    assert_equal(ret1.d, vals+vals)
    assert(type(ret1.units) is tuple)

    ret2 = index + row_index

    assert_equal(ret2.units, index.units)
    assert_equal(ret2.d, vals+row_vals)
    assert(type(ret2.units) is tuple)

    ret3 = row_index + index

    assert_equal(ret3.units, index.units)
    assert_equal(ret3.d, vals+row_vals)
    assert(type(ret3.units) is tuple)

def test_unit_conversions():
    vals = np.random.random((100, 3))
    index = IndexArray(vals, input_units=[u.km, u.g, u.s])

    converted = index.in_units(('m', 'kg', 'min'))
    correct_result = np.copy(vals)
    correct_result[:, 0] *= 1000
    correct_result[:, 1] /= 1000
    correct_result[:, 2] /= 60

    assert_equal(converted.units, (Unit('m'), Unit('kg'), Unit('min')))
    assert_almost_equal(correct_result, converted.d)
    assert(type(converted.units) is tuple)

    row_vals = np.random.random(3)
    row_index = IndexArray(row_vals, input_units=[u.km, u.g, u.s])

    converted = row_index.in_units(('m', 'kg', 'min'))
    correct_result = np.copy(row_vals)
    correct_result[0] *= 1000
    correct_result[1] /= 1000
    correct_result[2] /= 60

    assert_equal(converted.units, (Unit('m'), Unit('kg'), Unit('min')))
    assert_almost_equal(correct_result, converted.d)
    assert(type(converted.units) is tuple)

def compare_slicing(desired, actual, unit_type, array_type):
    assert_equal(desired, actual)
    assert_equal(desired.units, actual.units)
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
    ret7 = IndexArray(vals[[5, 7]], [u.km, u.g, u.s])

    sl1 = index[:, 0]
    sl2 = index[0, :]
    sl3 = index[..., 0]
    sl4 = index[1, 0]
    sl5 = index[12, 2]
    sl6 = index[5]
    sl7 = index[[5, 7]]

    compare_slicing(ret1, sl1, Unit, YTArray)
    compare_slicing(ret2, sl2, tuple, IndexArray)
    compare_slicing(ret3, sl3, Unit, YTArray)
    compare_slicing(ret4, sl4, Unit, YTQuantity)
    compare_slicing(ret5, sl5, Unit, YTQuantity)
    compare_slicing(ret6, sl6, tuple, IndexArray)
    compare_slicing(ret7, sl7, tuple, IndexArray)

    index = IndexArray([0, 1, 2], input_units=[u.km, u.g, u.s])
    compare_slicing(index[0], YTQuantity(0, 'km'), Unit, YTQuantity)
    compare_slicing(index[2], YTQuantity(2, 's'), Unit, YTQuantity)
    compare_slicing(index[:2], IndexArray([0, 1], [u.km, u.g]), tuple,
                    IndexArray)
    compare_slicing(index[1:2], YTQuantity(1, 'g'), Unit, YTQuantity)
    compare_slicing(index[:], index, tuple, IndexArray)
    compare_slicing(index[...], index, tuple, IndexArray)
    compare_slicing(index[[1, 2]], IndexArray([1, 2], [u.g, u.s]), tuple,
                    IndexArray)

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

def test_incompatible_unit_operations():
    index1 = IndexArray(np.random.random((6, 2)), [u.km, u.g])
    index2 = IndexArray(np.random.random((6, 3)), [u.km, u.s, u.g])

    assert_raises(ValueError, operator.add, index1, index2)
    assert_raises(ValueError, np.add, index1, index2)
