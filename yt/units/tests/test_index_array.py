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
from yt.units.index_array import YTIndexArray
from yt.units.yt_array import YTArray, YTQuantity
from yt.units.unit_object import Unit, UnitTuple

def test_creation():
    arr = YTArray([1, 2, 3], 'g')
    iarr = YTIndexArray(arr)

    assert_equal(iarr.units, (u.g.units, )*3)
    assert_equal(iarr.ndview, np.array([1, 2, 3]))
    assert(type(iarr.units) is UnitTuple)

    iarr = YTIndexArray([1, 2, 3], 'g')

    assert_equal(iarr.units, (u.g.units, )*3)
    assert_equal(iarr.ndview, np.array([1, 2, 3]))
    assert(type(iarr.units) is UnitTuple)

    iarr = YTIndexArray([1, 2, 3], Unit('g'))

    assert_equal(iarr.units, (u.g.units, )*3)
    assert_equal(iarr.ndview, np.array([1, 2, 3]))
    assert(type(iarr.units) is UnitTuple)

    iarr = YTIndexArray([1, 2, 3], ['g']*3)

    assert_equal(iarr.units, (u.g.units, )*3)
    assert_equal(iarr.ndview, np.array([1, 2, 3]))
    assert(type(iarr.units) is UnitTuple)

    iarr = YTIndexArray([1, 2, 3], [Unit('g')]*3)

    assert_equal(iarr.units, (u.g.units, )*3)
    assert_equal(iarr.ndview, np.array([1, 2, 3]))
    assert(type(iarr.units) is UnitTuple)

def test_creation_errors():
    assert_raises(RuntimeError, YTIndexArray, [1, 2, 3, 4], ['g', 'cm', 's'])
    assert_raises(NotImplementedError, YTIndexArray, np.random.random((10, 10, 3)),
                  'g')

def test_registry_association():
    ds = fake_random_ds(64)

    arr = YTIndexArray([1, 2, 3], ('g', 'km', 's'), registry=ds.unit_registry)

    assert(all([u.registry is ds.unit_registry for u in (3*arr).units]))
    assert(all([u.registry is ds.unit_registry for u in (arr*3).units]))

def assert_binary_op(oper1, oper2, operation, expected_value, expected_type,
                     expected_units):

    actual_result = operation(oper1, oper2)
    aunits = getattr(actual_result, 'units', None)

    assert_equal(actual_result.view(np.ndarray), expected_value)
    assert_equal(aunits, expected_units)
    assert(type(actual_result) is expected_type)
    assert(type(aunits) is type(expected_units))

def assert_commutative_binary_op(oper1, oper2, operation, expected_value,
                                 expected_type, expected_units):

    assert_binary_op(oper1, oper2, operation, expected_value, expected_type,
                     expected_units)
    assert_binary_op(oper2, oper1, operation, expected_value, expected_type,
                     expected_units)

def test_multiplication():
    vals = np.random.random((100, 3))
    index = YTIndexArray(vals, input_units=[u.km, u.g, u.s])
    row_index = index[0]

    for operand, value in zip((index, row_index), (vals, vals[0])):
        for operation in (np.multiply, operator.mul):
            assert_commutative_binary_op(
                operand, index, operation, value * vals, YTIndexArray,
                UnitTuple((u.km**2).units, (u.g**2).units, (u.s**2).units))
            assert_commutative_binary_op(
                operand, u.km, operation, value, YTIndexArray,
                UnitTuple((u.km**2).units, (u.g*u.km).units, (u.s*u.km).units))
            assert_commutative_binary_op(
                operand, 2, operation, 2 * value, YTIndexArray, index.units)

def test_division():
    vals = np.random.random((100, 3))
    index = YTIndexArray(vals, input_units=[u.km, u.g, u.s])
    row_index = index[0]

    for operand, value in zip((index, row_index), (vals, vals[0])):
        for operation in (np.divide, operator.div):
            assert_binary_op(
                operand, index, operation, value / vals, YTIndexArray,
                UnitTuple(Unit(), Unit(), Unit()))
            assert_binary_op(
                index, operand, operation, vals / value, YTIndexArray,
                UnitTuple(Unit(), Unit(), Unit()))

            assert_binary_op(
                operand, u.km, operation, value, YTIndexArray,
                UnitTuple(Unit(), (u.g/u.km).units, (u.s/u.km).units))
            assert_binary_op(
                u.km, operand, operation, 1/value, YTIndexArray,
                UnitTuple(Unit(), (u.km/u.g).units, (u.km/u.s).units))

            assert_binary_op(
                operand, 2, operation, value / 2, YTIndexArray, index.units)
            assert_binary_op(
                2, operand, operation, 2 / value, YTIndexArray,
                UnitTuple(Unit('1/km'), Unit('1/g'), Unit('1/s')))

def test_addition():
    vals = np.random.random((100, 3))
    index = YTIndexArray(vals, input_units=[u.km, u.g, u.s])

    row_vals = np.random.random(3)
    row_index = YTIndexArray(row_vals, input_units=[u.km, u.g, u.s])

    for operand, value in zip((index, row_index), (vals, row_vals)):
        for operation in (np.add, operator.add):
            assert_commutative_binary_op(
                operand, index, operation, value + vals, YTIndexArray, index.units)
            assert_commutative_binary_op(
                operand, row_index, operation, value + row_vals, YTIndexArray,
                index.units)

def test_subtraction():
    vals = np.random.random((100, 3))
    index = YTIndexArray(vals, input_units=[u.km, u.g, u.s])

    row_vals = np.random.random(3)
    row_index = YTIndexArray(row_vals, input_units=[u.km, u.g, u.s])

    for operand, value in zip((index, row_index), (vals, row_vals)):
        for operation in (np.subtract, operator.sub):
            assert_binary_op(
                operand, index, operation, value - vals, YTIndexArray, index.units)
            assert_binary_op(
                index, operand, operation, vals - value, YTIndexArray, index.units)
            assert_binary_op(
                operand, row_index, operation, value - row_vals, YTIndexArray,
                index.units)
            assert_binary_op(
                row_index, operand, operation, row_vals - value, YTIndexArray,
                index.units)

def test_comparisons():
    vals = np.random.random((100, 3))
    index = YTIndexArray(vals, input_units=[u.km, u.km, u.km])
    arr = YTArray(vals, input_units=u.km)

    for operation in (np.equal, operator.eq, np.greater, operator.gt,
                      np.less, operator.lt, np.greater_equal, operator.ge,
                      np.less_equal, operator.le):
        assert_commutative_binary_op(
            index, arr, operation, operation(vals, vals), np.ndarray, None)

def test_homogenous_unit_operations():
    vals = np.random.random((100, 3))
    index = YTIndexArray(vals, input_units=[u.km, u.km, u.km])
    row_index = index[0]
    arr = YTArray(vals, input_units=u.km)
    row_arr = arr[0]
    quantity = row_arr[0]

    assert_equal(index, arr)

    for operand, value in zip((index, row_index), (vals, vals[0])):
        for operation in (np.add, operator.add):
            assert_commutative_binary_op(
                operand, arr, operation, value + vals, YTIndexArray, index.units)
            assert_commutative_binary_op(
                operand, row_arr, operation, value + vals[0], YTIndexArray,
                index.units)
            assert_commutative_binary_op(
                operand, quantity, operation, value + vals[0, 0], YTIndexArray,
                index.units)
        for operation in (np.multiply, operator.mul):
            assert_commutative_binary_op(
                operand, arr, operation, value * vals, YTIndexArray,
                UnitTuple(((u.km**2).units, )*3))
            assert_commutative_binary_op(
                operand, row_arr, operation, value * vals[0], YTIndexArray,
                UnitTuple(((u.km**2).units, )*3))
            assert_commutative_binary_op(
                operand, quantity, operation, value * vals[0, 0], YTIndexArray,
                UnitTuple(((u.km**2).units, )*3))

def test_unit_conversions():
    vals = np.random.random((100, 3))
    index = YTIndexArray(vals, input_units=[u.km, u.g, u.s])

    converted = index.in_units(('m', 'kg', 'min'))
    correct_result = np.copy(vals)
    correct_result[:, 0] *= 1000
    correct_result[:, 1] /= 1000
    correct_result[:, 2] /= 60

    assert_equal(converted.units, UnitTuple(Unit('m'), Unit('kg'), Unit('min')))
    assert_almost_equal(correct_result, converted.d)
    assert(type(converted.units) is UnitTuple)

    row_vals = np.random.random(3)
    row_index = YTIndexArray(row_vals, input_units=[u.km, u.g, u.s])

    converted = row_index.in_units(('m', 'kg', 'min'))
    correct_result = np.copy(row_vals)
    correct_result[0] *= 1000
    correct_result[1] /= 1000
    correct_result[2] /= 60

    assert_equal(converted.units, UnitTuple(Unit('m'), Unit('kg'), Unit('min')))
    assert_almost_equal(correct_result, converted.d)
    assert(type(converted.units) is UnitTuple)

def compare_slicing(desired, actual, unit_type, array_type):
    assert_equal(desired, actual)
    assert_equal(desired.units, actual.units)
    assert(type(actual.units) is unit_type)
    assert(type(actual) is array_type)

def test_slicing():
    vals = np.arange(300)
    vals.shape = (100, 3)
    index = YTIndexArray(vals, input_units=[u.km, u.g, u.s])

    ret1 = YTArray(3*np.arange(100), 'km')
    ret2 = YTIndexArray(np.arange(3), [u.km, u.g, u.s])
    ret3 = ret1
    ret4 = YTQuantity(3, 'km')
    ret5 = YTQuantity(38, 's')
    ret6 = YTIndexArray([15, 16, 17], [u.km, u.g, u.s])
    ret7 = YTIndexArray(vals[[5, 7]], [u.km, u.g, u.s])

    sl1 = index[:, 0]
    sl2 = index[0, :]
    sl3 = index[..., 0]
    sl4 = index[1, 0]
    sl5 = index[12, 2]
    sl6 = index[5]
    sl7 = index[[5, 7]]

    compare_slicing(ret1, sl1, Unit, YTArray)
    compare_slicing(ret2, sl2, UnitTuple, YTIndexArray)
    compare_slicing(ret3, sl3, Unit, YTArray)
    compare_slicing(ret4, sl4, Unit, YTQuantity)
    compare_slicing(ret5, sl5, Unit, YTQuantity)
    compare_slicing(ret6, sl6, UnitTuple, YTIndexArray)
    compare_slicing(ret7, sl7, UnitTuple, YTIndexArray)

    index = YTIndexArray([0, 1, 2], input_units=[u.km, u.g, u.s])
    compare_slicing(index[0], YTQuantity(0, 'km'), Unit, YTQuantity)
    compare_slicing(index[2], YTQuantity(2, 's'), Unit, YTQuantity)
    compare_slicing(index[:2], YTIndexArray([0, 1], [u.km, u.g]), UnitTuple,
                    YTIndexArray)
    compare_slicing(index[1:2], YTQuantity(1, 'g'), Unit, YTQuantity)
    compare_slicing(index[:], index, UnitTuple, YTIndexArray)
    compare_slicing(index[...], index, UnitTuple, YTIndexArray)
    compare_slicing(index[[1, 2]], YTIndexArray([1, 2], [u.g, u.s]), UnitTuple,
                    YTIndexArray)

def test_boolean_unary_ops():
    vals = np.arange(30)
    vals.shape = (10, 3)
    index = YTIndexArray(vals, [u.km, u.km, u.km])

    assert_equal(np.isnan(index), np.isnan(vals))
    assert(all([i is np.bool_(False) for i in np.isnan(index).ravel()]))

def test_str_and_repr():
    vals = np.arange(6)
    vals.shape = (2, 3)
    index = YTIndexArray(vals, input_units=[u.km, u.g, u.s])

    assert_equal(str(index), '[[ 0.  1.  2.]\n [ 3.  4.  5.]] (km, g, s)')
    assert_equal(
        repr(index),
        'YTIndexArray([[ 0.,  1.,  2.],\n       [ 3.,  4.,  5.]]) (km, g, s)')

def test_incompatible_unit_operations():
    index1 = YTIndexArray(np.random.random((6, 2)), [u.km, u.g])
    index2 = YTIndexArray(np.random.random((6, 3)), [u.km, u.s, u.g])

    assert_raises(ValueError, operator.add, index1, index2)
    assert_raises(ValueError, np.add, index1, index2)
