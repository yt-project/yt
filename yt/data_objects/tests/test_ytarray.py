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
from yt.data_objects.yt_array import YTArray, YTQuantity
from yt.utilities.exceptions import YTUnitOperationError
from yt.funcs import fix_length
import operator

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

    # different units
    a1 = YTArray([1, 2, 3], 'cm')
    a2 = YTArray([4, 5, 6], 'm')
    answer1 = YTArray([401, 502, 603], 'cm')
    answer2 = YTArray([4.01, 5.02, 6.03], 'm')
    
    yield operate_and_compare, a1, a2, operator.add, answer1
    yield operate_and_compare, a2, a1, operator.add, answer2

    # Test dimensionless quantities
    a1 = YTArray([1,2,3])
    a2 = array([4,5,6])
    answer = YTArray([5, 7, 9])

    yield operate_and_compare, a1, a2, operator.add, answer
    yield operate_and_compare, a2, a1, operator.add, answer

    # Catch the different dimensions error
    a1 = YTArray([1, 2, 3], 'm')
    a2 = YTArray([4, 5, 6], 'kg')

    try:
        a1+a2
    except YTUnitOperationError as e:
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

    # different units
    a1 = YTArray([1, 2, 3], 'cm')
    a2 = YTArray([4, 5, 6], 'm')
    answer1 = YTArray([-399, -498, -597], 'cm')
    answer2 = YTArray([3.99, 4.98, 5.97], 'm')

    yield operate_and_compare, a1, a2, operator.sub, answer1
    yield operate_and_compare, a2, a1, operator.sub, answer2

    # Test dimensionless quantities
    a1 = YTArray([1,2,3])
    a2 = array([4,5,6])
    answer1 = YTArray([-3, -3, -3])
    answer2 = YTArray([3, 3, 3])

    yield operate_and_compare, a1, a2, operator.sub, answer1
    yield operate_and_compare, a2, a1, operator.sub, answer2

    # Catch the different dimensions error
    a1 = YTArray([1, 2, 3], 'm')
    a2 = YTArray([4, 5, 6], 'kg')

    try:
        a1-a2
    except YTUnitOperationError as e:
        assert isinstance(e, YTUnitOperationError)
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

    # different units, same dimension
    a1 = YTArray([1, 2, 3], 'cm')
    a2 = YTArray([4, 5, 6], 'm')
    answer1 = YTArray([400, 1000, 1800], 'cm**2')
    answer2 = YTArray([.04, .10, .18], 'm**2')
    
    yield operate_and_compare, a1, a2, operator.mul, answer1
    yield operate_and_compare, a2, a1, operator.mul, answer2

    # different dimensions
    a1 = YTArray([1, 2, 3], 'cm')
    a2 = YTArray([4, 5, 6], 'g')
    answer = YTArray([4, 10, 18], 'cm*g')
        
    yield operate_and_compare, a1, a2, operator.mul, answer
    yield operate_and_compare, a2, a1, operator.mul, answer

    # One dimensionless, one unitful
    a1 = YTArray([1,2,3], 'cm')
    a2 = array([4,5,6])
    answer = YTArray([4, 10, 18], 'cm')

    yield operate_and_compare, a1, a2, operator.mul, answer
    yield operate_and_compare, a2, a1, operator.mul, answer

    # Both dimensionless quantities
    a1 = YTArray([1,2,3])
    a2 = array([4,5,6])
    answer = YTArray([4, 10, 18])

    yield operate_and_compare, a1, a2, operator.mul, answer
    yield operate_and_compare, a2, a1, operator.mul, answer

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

    # different units, same dimension
    a1 = YTArray([1., 2., 3.], 'cm')
    a2 = YTArray([4., 5., 6.], 'm')
    answer1 = YTArray([.0025, .004, .005])
    answer2 = YTArray([400, 250, 200])
    
    yield operate_and_compare, a1, a2, operator.div, answer1
    yield operate_and_compare, a2, a1, operator.div, answer2

    # different dimensions
    a1 = YTArray([1., 2., 3.], 'cm')
    a2 = YTArray([4., 5., 6.], 'g')
    answer1 = YTArray([0.25, 0.4, 0.5], 'cm/g')
    answer2 = YTArray([4, 2.5, 2], 'g/cm')

    yield operate_and_compare, a1, a2, operator.div, answer1
    yield operate_and_compare, a2, a1, operator.div, answer2

    # One dimensionless, one unitful
    a1 = YTArray([1., 2., 3.], 'cm')
    a2 = array([4., 5., 6.])
    answer1 = YTArray([0.25, 0.4, 0.5], 'cm')
    answer2 = YTArray([4, 2.5, 2], '1/cm')

    yield operate_and_compare, a1, a2, operator.div, answer1
    yield operate_and_compare, a2, a1, operator.div, answer2

    # Both dimensionless quantities
    a1 = YTArray([1., 2., 3.])
    a2 = array([4., 5., 6.])
    answer1 = YTArray([0.25, 0.4, 0.5])
    answer2 = YTArray([4, 2.5, 2])

    yield operate_and_compare, a1, a2, operator.div, answer1
    yield operate_and_compare, a2, a1, operator.div, answer2

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
    length = YTQuantity(1.0,'code_length')
    new_length = fix_length(length)
    assert length == new_length
