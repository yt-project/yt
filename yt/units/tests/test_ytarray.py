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

from numpy.testing import assert_approx_equal, assert_array_equal
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
    try:
        np.add(a1, a2)
    except YTUfuncUnitError:
        pass

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

    try:
        a1+a2
    except YTUnitOperationError as e:
        pass

    assert(isinstance(e, YTUnitOperationError))
    assert str(e) == \
        'The addition operator for YTArrays with units ' \
        '(m) and (kg) is not well defined.'

    try:
        a1 += a2
    except YTUnitOperationError as e:
        pass

    assert(isinstance(e, YTUnitOperationError))
    assert str(e) == \
        'The addition operator for YTArrays with units ' \
        '(m) and (kg) is not well defined.'

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

    try:
        np.subtract(a1, a2)
    except YTUfuncUnitError:
        pass

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

    try:
        a1-a2
    except YTUnitOperationError as e:
        pass

    assert isinstance(e, YTUnitOperationError)
    assert str(e) == \
        'The subtraction operator for YTArrays with units ' \
        '(m) and (kg) is not well defined.'

    try:
        a1 -= a2
    except YTUnitOperationError as e:
        pass

    assert(isinstance(e, YTUnitOperationError))
    assert str(e) == \
        'The subtraction operator for YTArrays with units ' \
        '(m) and (kg) is not well defined.'

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

def test_yt_array_yt_quantity_ops():
    """
    Test operations that combine YTArray and YTQuantity
    """
    a = YTArray(range(10), 'cm')
    b = YTQuantity(5, 'g')

    assert isinstance(a*b, YTArray)
    assert isinstance(b*a, YTArray)

    assert isinstance(a/b, YTArray)
    assert isinstance(b/a, YTArray)

    assert isinstance(a*a, YTArray)
    assert isinstance(a/a, YTArray)

    assert isinstance(b*b, YTQuantity)
    assert isinstance(b/b, YTQuantity)

def test_selecting():
    """
    Test slicing of two YTArrays

    """
    a = YTArray(range(10), 'cm')
    assert_array_equal, a[:3], YTArray([0, 1, 2], 'cm')
    assert isinstance(a[0], YTQuantity)

def test_fix_length():
    """
    Test fixing the length of an array. Used in spheres and other data objects
    """
    pf = fake_random_pf(64, nprocs=1)
    length = pf.quan(1.0,'code_length')
    new_length = fix_length(length, pf=pf)
    assert length == new_length

def test_ytarray_pickle():
    pf = fake_random_pf(64, nprocs=1)
    test_data = [pf.quan(12.0, 'code_length'), pf.arr([1,2,3], 'code_length')]
    
    for data in test_data:
        tempf = tempfile.NamedTemporaryFile(delete=False)
        pickle.dump(data, tempf)
        tempf.close()

        loaded_data = pickle.load(open(tempf.name, "rb"))

        assert_array_equal(data, loaded_data)
        assert data.units == loaded_data.units
        assert_array_equal(array(data.in_cgs()), array(loaded_data.in_cgs()))
        assert float(data.units.cgs_value) == float(loaded_data.units.cgs_value)

def test_deepcopy():
    quan = YTQuantity(1, 'g')
    arr = YTArray([1,2,3], 'cm')

    assert copy.deepcopy(quan) == quan
    assert_array_equal(copy.deepcopy(arr), arr)
