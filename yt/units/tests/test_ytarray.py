"""
Test ndarray subclass that handles symbolic units.




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from nose.tools import assert_true
from numpy.testing import \
    assert_approx_equal, assert_array_equal, \
    assert_equal, assert_raises
from numpy import array
from yt.units.yt_array import YTArray, YTQuantity
from yt.utilities.exceptions import YTUnitOperationError, YTUfuncUnitError
from yt.testing import fake_random_pf
from yt.funcs import fix_length
import numpy as np
import copy
import operator
import cPickle as pickle
import tempfile

def operate_and_compare(a, b, op, answer):
    # Test generator for YTArrays tests
    assert_array_equal(op(a, b), answer)

def assert_isinstance(a, type):
    assert isinstance(a, type)

def test_addition():
    """
    Test addition of two YTArrays

    """

    # Same units
    a1 = YTArray([1, 2, 3], 'cm')
    a2 = YTArray([4, 5, 6], 'cm')
    answer = YTArray([5, 7, 9], 'cm')

    yield operate_and_compare, a1, a2, operator.add, answer
    yield operate_and_compare, a2, a1, operator.add, answer
    yield operate_and_compare, a1, a2, np.add, answer
    yield operate_and_compare, a2, a1, np.add, answer

    # different units
    a1 = YTArray([1, 2, 3], 'cm')
    a2 = YTArray([4, 5, 6], 'm')
    answer1 = YTArray([401, 502, 603], 'cm')
    answer2 = YTArray([4.01, 5.02, 6.03], 'm')

    yield operate_and_compare, a1, a2, operator.add, answer1
    yield operate_and_compare, a2, a1, operator.add, answer2
    yield assert_raises, YTUfuncUnitError, np.add, a1, a2

    # Test dimensionless quantities
    a1 = YTArray([1,2,3])
    a2 = array([4,5,6])
    answer = YTArray([5, 7, 9])

    yield operate_and_compare, a1, a2, operator.add, answer
    yield operate_and_compare, a2, a1, operator.add, answer
    yield operate_and_compare, a1, a2, np.add, answer
    yield operate_and_compare, a2, a1, np.add, answer

    # Catch the different dimensions error
    a1 = YTArray([1, 2, 3], 'm')
    a2 = YTArray([4, 5, 6], 'kg')

    yield assert_raises, YTUnitOperationError, operator.add, a1, a2
    yield assert_raises, YTUnitOperationError, operator.iadd, a1, a2

def test_subtraction():
    """
    Test subtraction of two YTArrays

    """

    # Same units
    a1 = YTArray([1, 2, 3], 'cm')
    a2 = YTArray([4, 5, 6], 'cm')
    answer1 = YTArray([-3, -3, -3], 'cm')
    answer2 = YTArray([3, 3, 3], 'cm')

    yield operate_and_compare, a1, a2, operator.sub, answer1
    yield operate_and_compare, a2, a1, operator.sub, answer2
    yield operate_and_compare, a1, a2, np.subtract, answer1
    yield operate_and_compare, a2, a1, np.subtract, answer2

    # different units
    a1 = YTArray([1, 2, 3], 'cm')
    a2 = YTArray([4, 5, 6], 'm')
    answer1 = YTArray([-399, -498, -597], 'cm')
    answer2 = YTArray([3.99, 4.98, 5.97], 'm')

    yield operate_and_compare, a1, a2, operator.sub, answer1
    yield operate_and_compare, a2, a1, operator.sub, answer2
    yield assert_raises, YTUfuncUnitError, np.subtract, a1, a2

    # Test dimensionless quantities
    a1 = YTArray([1,2,3])
    a2 = array([4,5,6])
    answer1 = YTArray([-3, -3, -3])
    answer2 = YTArray([3, 3, 3])

    yield operate_and_compare, a1, a2, operator.sub, answer1
    yield operate_and_compare, a2, a1, operator.sub, answer2
    yield operate_and_compare, a1, a2, np.subtract, answer1
    yield operate_and_compare, a2, a1, np.subtract, answer2

    # Catch the different dimensions error
    a1 = YTArray([1, 2, 3], 'm')
    a2 = YTArray([4, 5, 6], 'kg')

    yield assert_raises, YTUnitOperationError, operator.sub, a1, a2
    yield assert_raises, YTUnitOperationError, operator.isub, a1, a2

def test_multiplication():
    """
    Test multiplication of two YTArrays

    """

    # Same units
    a1 = YTArray([1, 2, 3], 'cm')
    a2 = YTArray([4, 5, 6], 'cm')
    answer = YTArray([4, 10, 18], 'cm**2')

    yield operate_and_compare, a1, a2, operator.mul, answer
    yield operate_and_compare, a2, a1, operator.mul, answer
    yield operate_and_compare, a1, a2, np.multiply, answer
    yield operate_and_compare, a2, a1, np.multiply, answer

    # different units, same dimension
    a1 = YTArray([1, 2, 3], 'cm')
    a2 = YTArray([4, 5, 6], 'm')
    answer1 = YTArray([400, 1000, 1800], 'cm**2')
    answer2 = YTArray([.04, .10, .18], 'm**2')
    answer3 = YTArray([4, 10, 18], 'cm*m')

    yield operate_and_compare, a1, a2, operator.mul, answer1
    yield operate_and_compare, a2, a1, operator.mul, answer2
    yield operate_and_compare, a1, a2, np.multiply, answer3
    yield operate_and_compare, a2, a1, np.multiply, answer3

    # different dimensions
    a1 = YTArray([1, 2, 3], 'cm')
    a2 = YTArray([4, 5, 6], 'g')
    answer = YTArray([4, 10, 18], 'cm*g')

    yield operate_and_compare, a1, a2, operator.mul, answer
    yield operate_and_compare, a2, a1, operator.mul, answer
    yield operate_and_compare, a1, a2, np.multiply, answer
    yield operate_and_compare, a2, a1, np.multiply, answer

    # One dimensionless, one unitful
    a1 = YTArray([1,2,3], 'cm')
    a2 = array([4,5,6])
    answer = YTArray([4, 10, 18], 'cm')

    yield operate_and_compare, a1, a2, operator.mul, answer
    yield operate_and_compare, a2, a1, operator.mul, answer
    yield operate_and_compare, a1, a2, np.multiply, answer
    yield operate_and_compare, a2, a1, np.multiply, answer

    # Both dimensionless quantities
    a1 = YTArray([1,2,3])
    a2 = array([4,5,6])
    answer = YTArray([4, 10, 18])

    yield operate_and_compare, a1, a2, operator.mul, answer
    yield operate_and_compare, a2, a1, operator.mul, answer
    yield operate_and_compare, a1, a2, np.multiply, answer
    yield operate_and_compare, a2, a1, np.multiply, answer

def test_division():
    """
    Test multiplication of two YTArrays

    """

    # Same units
    a1 = YTArray([1., 2., 3.], 'cm')
    a2 = YTArray([4., 5., 6.], 'cm')
    answer1 = YTArray([0.25, 0.4, 0.5])
    answer2 = YTArray([4, 2.5, 2])

    yield operate_and_compare, a1, a2, operator.div, answer1
    yield operate_and_compare, a2, a1, operator.div, answer2
    yield operate_and_compare, a1, a2, np.divide, answer1
    yield operate_and_compare, a2, a1, np.divide, answer2

    # different units, same dimension
    a1 = YTArray([1., 2., 3.], 'cm')
    a2 = YTArray([4., 5., 6.], 'm')
    answer1 = YTArray([.0025, .004, .005])
    answer2 = YTArray([400, 250, 200])
    answer3 = YTArray([0.25, 0.4, 0.5], 'cm/m')
    answer4 = YTArray([4.0, 2.5, 2.0], 'm/cm')

    yield operate_and_compare, a1, a2, operator.div, answer1
    yield operate_and_compare, a2, a1, operator.div, answer2
    yield operate_and_compare, a1, a2, np.divide, answer3
    yield operate_and_compare, a2, a1, np.divide, answer4

    # different dimensions
    a1 = YTArray([1., 2., 3.], 'cm')
    a2 = YTArray([4., 5., 6.], 'g')
    answer1 = YTArray([0.25, 0.4, 0.5], 'cm/g')
    answer2 = YTArray([4, 2.5, 2], 'g/cm')

    yield operate_and_compare, a1, a2, operator.div, answer1
    yield operate_and_compare, a2, a1, operator.div, answer2
    yield operate_and_compare, a1, a2, np.divide, answer1
    yield operate_and_compare, a2, a1, np.divide, answer2

    # One dimensionless, one unitful
    a1 = YTArray([1., 2., 3.], 'cm')
    a2 = array([4., 5., 6.])
    answer1 = YTArray([0.25, 0.4, 0.5], 'cm')
    answer2 = YTArray([4, 2.5, 2], '1/cm')

    yield operate_and_compare, a1, a2, operator.div, answer1
    yield operate_and_compare, a2, a1, operator.div, answer2
    yield operate_and_compare, a1, a2, np.divide, answer1
    yield operate_and_compare, a2, a1, np.divide, answer2

    # Both dimensionless quantities
    a1 = YTArray([1., 2., 3.])
    a2 = array([4., 5., 6.])
    answer1 = YTArray([0.25, 0.4, 0.5])
    answer2 = YTArray([4, 2.5, 2])

    yield operate_and_compare, a1, a2, operator.div, answer1
    yield operate_and_compare, a2, a1, operator.div, answer2
    yield operate_and_compare, a1, a2, np.divide, answer1
    yield operate_and_compare, a2, a1, np.divide, answer2

def test_comparisons():
    """
    Test numpy ufunc comparison operators for unit consistency.

    """
    from yt.units.yt_array import YTArray, YTQuantity

    a1 = YTArray([1,2,3], 'cm')
    a2 = YTArray([2,1,3], 'cm')
    a3 = YTArray([.02, .01, .03], 'm')

    ops = (
        np.less,
        np.less_equal,
        np.greater,
        np.greater_equal,
        np.equal,
        np.not_equal
    )

    answers = (
        [True, False, False],
        [True, False, True],
        [False, True, False],
        [False, True, True],
        [False, False, True],
        [True, True, False],
    )

    for op, answer in zip(ops, answers):
        yield operate_and_compare, a1, a2, op, answer

    for op in ops:
        yield assert_raises, YTUfuncUnitError, op, a1, a3

    for op, answer in zip(ops, answers):
        yield operate_and_compare, a1, a3.in_units('cm'), op, answer

def test_unit_conversions():
    """
    Test operations that convert to different units or cast to ndarray

    """
    from yt.units.yt_array import YTQuantity
    from yt.units.unit_object import Unit

    km = YTQuantity(1, 'km')
    km_in_cm = km.in_units('cm')
    km_unit = Unit('km')
    cm_unit = Unit('cm')

    yield assert_equal, km_in_cm, km
    yield assert_equal, km_in_cm.in_cgs(), 1e5
    yield assert_equal, km_in_cm.units, cm_unit

    km.convert_to_units('cm')

    yield assert_equal, km, YTQuantity(1, 'km')
    yield assert_equal, km.in_cgs(), 1e5
    yield assert_equal, km.units, cm_unit

    yield assert_isinstance, km.to_ndarray(), np.ndarray
    yield assert_isinstance, km.ndarray_view(), np.ndarray


def test_yt_array_yt_quantity_ops():
    """
    Test operations that combine YTArray and YTQuantity
    """
    a = YTArray(range(10), 'cm')
    b = YTQuantity(5, 'g')

    yield assert_isinstance, a*b, YTArray
    yield assert_isinstance, b*a, YTArray

    yield assert_isinstance, a/b, YTArray
    yield assert_isinstance, b/a, YTArray

    yield assert_isinstance, a*a, YTArray
    yield assert_isinstance, a/a, YTArray

    yield assert_isinstance, b*b, YTQuantity
    yield assert_isinstance, b/b, YTQuantity

def test_selecting():
    """
    Test slicing of two YTArrays

    """
    a = YTArray(range(10), 'cm')
    a_slice = a[:3]
    a_fancy_index = a[[1,1,3,5]]
    a_array_fancy_index = a[array([[1,1], [3,5]])]
    a_boolean_index = a[a > 5]
    a_selection = a[0]

    yield assert_array_equal, a_slice, YTArray([0, 1, 2], 'cm')
    yield assert_array_equal, a_fancy_index, YTArray([1,1,3,5], 'cm')
    yield assert_array_equal, a_array_fancy_index, YTArray([[1, 1,], [3,5]], 'cm')
    yield assert_array_equal, a_boolean_index, YTArray([6,7,8,9], 'cm')
    yield assert_isinstance, a_selection, YTQuantity

    # .base points to the original array for a numpy view.  If it is not a
    # view, .base is None.
    yield assert_true, a_slice.base is a
    yield assert_true, a_fancy_index.base is None
    yield assert_true, a_array_fancy_index.base is None
    yield assert_true, a_boolean_index.base is None

def test_fix_length():
    """
    Test fixing the length of an array. Used in spheres and other data objects
    """
    pf = fake_random_pf(64, nprocs=1)
    length = pf.quan(1.0,'code_length')
    new_length = fix_length(length, pf=pf)
    yield assert_equal, length, new_length

def test_ytarray_pickle():
    pf = fake_random_pf(64, nprocs=1)
    test_data = [pf.quan(12.0, 'code_length'), pf.arr([1,2,3], 'code_length')]
    
    for data in test_data:
        tempf = tempfile.NamedTemporaryFile(delete=False)
        pickle.dump(data, tempf)
        tempf.close()

        loaded_data = pickle.load(open(tempf.name, "rb"))

        yield assert_array_equal, data, loaded_data
        yield assert_equal, data.units, loaded_data.units
        yield assert_array_equal, array(data.in_cgs()), array(loaded_data.in_cgs())
        yield assert_equal, float(data.units.cgs_value), float(loaded_data.units.cgs_value)

def test_copy():
    quan = YTQuantity(1, 'g')
    arr = YTArray([1,2,3], 'cm')

    yield assert_equal, copy.copy(quan), quan
    yield assert_array_equal, copy.copy(arr), arr

    yield assert_equal,  copy.deepcopy(quan), quan
    yield assert_array_equal, copy.deepcopy(arr), arr

    yield assert_equal, quan.copy(), quan
    yield assert_array_equal, arr.copy(), arr

    yield assert_equal, np.copy(quan), quan
    yield assert_array_equal, np.copy(arr), arr
