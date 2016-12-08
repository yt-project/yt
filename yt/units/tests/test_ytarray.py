"""
Test ndarray subclass that handles symbolic units.




"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import copy
from yt.extern.six.moves import cPickle as pickle
import itertools
import numpy as np
import operator
import os
import shutil
import tempfile

from nose.tools import assert_true
from numpy.testing import \
    assert_array_equal, \
    assert_equal, assert_raises, \
    assert_array_almost_equal_nulp, \
    assert_array_almost_equal, \
    assert_almost_equal
from numpy import array
from yt.units.yt_array import \
    YTArray, YTQuantity, \
    unary_operators, binary_operators, \
    uconcatenate, uintersect1d, \
    uunion1d, loadtxt, savetxt
from yt.utilities.exceptions import \
    YTUnitOperationError, YTUfuncUnitError
from yt.testing import \
    fake_random_ds, requires_module, \
    assert_allclose_units
from yt.funcs import fix_length
from yt.units.unit_symbols import \
    cm, m, g, degree
from yt.utilities.physical_ratios import \
    metallicity_sun

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
    a3 = [4*cm, 5*cm, 6*cm]
    answer = YTArray([5, 7, 9], 'cm')

    yield operate_and_compare, a1, a2, operator.add, answer
    yield operate_and_compare, a2, a1, operator.add, answer
    yield operate_and_compare, a1, a3, operator.add, answer
    yield operate_and_compare, a3, a1, operator.add, answer
    yield operate_and_compare, a2, a1, np.add, answer
    yield operate_and_compare, a1, a2, np.add, answer
    yield operate_and_compare, a1, a3, np.add, answer
    yield operate_and_compare, a3, a1, np.add, answer

    # different units
    a1 = YTArray([1, 2, 3], 'cm')
    a2 = YTArray([4, 5, 6], 'm')
    a3 = [4*m, 5*m, 6*m]
    answer1 = YTArray([401, 502, 603], 'cm')
    answer2 = YTArray([4.01, 5.02, 6.03], 'm')

    yield operate_and_compare, a1, a2, operator.add, answer1
    yield operate_and_compare, a2, a1, operator.add, answer2
    yield operate_and_compare, a1, a3, operator.add, answer1
    yield operate_and_compare, a3, a1, operator.add, answer1
    yield assert_raises, YTUfuncUnitError, np.add, a1, a2
    yield assert_raises, YTUfuncUnitError, np.add, a1, a3

    # Test dimensionless quantities
    a1 = YTArray([1, 2, 3])
    a2 = array([4, 5, 6])
    a3 = [4, 5, 6]
    answer = YTArray([5, 7, 9])

    yield operate_and_compare, a1, a2, operator.add, answer
    yield operate_and_compare, a2, a1, operator.add, answer
    yield operate_and_compare, a1, a3, operator.add, answer
    yield operate_and_compare, a3, a1, operator.add, answer
    yield operate_and_compare, a1, a2, np.add, answer
    yield operate_and_compare, a2, a1, np.add, answer
    yield operate_and_compare, a1, a3, np.add, answer
    yield operate_and_compare, a3, a1, np.add, answer

    # Catch the different dimensions error
    a1 = YTArray([1, 2, 3], 'm')
    a2 = YTArray([4, 5, 6], 'kg')

    yield assert_raises, YTUnitOperationError, operator.add, a1, a2
    yield assert_raises, YTUnitOperationError, operator.iadd, a1, a2

    # adding with zero is allowed irrespective of the units
    zeros = np.zeros(3)
    zeros_yta_dimless = YTArray(zeros, 'dimensionless')
    zeros_yta_length = YTArray(zeros, 'm')
    zeros_yta_mass = YTArray(zeros, 'kg')
    operands = [0, YTQuantity(0), YTQuantity(0, 'kg'), zeros, zeros_yta_dimless,
                zeros_yta_length, zeros_yta_mass]

    for op in [operator.add, np.add]:
        for operand in operands:
            yield operate_and_compare, a1, operand, op, a1
            yield operate_and_compare, operand, a1, op, a1
            yield operate_and_compare, 4*m, operand, op, 4*m
            yield operate_and_compare, operand, 4*m, op, 4*m

def test_subtraction():
    """
    Test subtraction of two YTArrays

    """

    # Same units
    a1 = YTArray([1, 2, 3], 'cm')
    a2 = YTArray([4, 5, 6], 'cm')
    a3 = [4*cm, 5*cm, 6*cm]
    answer1 = YTArray([-3, -3, -3], 'cm')
    answer2 = YTArray([3, 3, 3], 'cm')

    yield operate_and_compare, a1, a2, operator.sub, answer1
    yield operate_and_compare, a2, a1, operator.sub, answer2
    yield operate_and_compare, a1, a3, operator.sub, answer1
    yield operate_and_compare, a3, a1, operator.sub, answer2
    yield operate_and_compare, a1, a2, np.subtract, answer1
    yield operate_and_compare, a2, a1, np.subtract, answer2
    yield operate_and_compare, a1, a3, np.subtract, answer1
    yield operate_and_compare, a3, a1, np.subtract, answer2

    # different units
    a1 = YTArray([1, 2, 3], 'cm')
    a2 = YTArray([4, 5, 6], 'm')
    a3 = [4*m, 5*m, 6*m]
    answer1 = YTArray([-399, -498, -597], 'cm')
    answer2 = YTArray([3.99, 4.98, 5.97], 'm')
    answer3 = YTArray([399, 498, 597], 'cm')

    yield operate_and_compare, a1, a2, operator.sub, answer1
    yield operate_and_compare, a2, a1, operator.sub, answer2
    yield operate_and_compare, a1, a3, operator.sub, answer1
    yield operate_and_compare, a3, a1, operator.sub, answer3
    yield assert_raises, YTUfuncUnitError, np.subtract, a1, a2
    yield assert_raises, YTUfuncUnitError, np.subtract, a1, a3

    # Test dimensionless quantities
    a1 = YTArray([1, 2, 3])
    a2 = array([4, 5, 6])
    a3 = [4, 5, 6]
    answer1 = YTArray([-3, -3, -3])
    answer2 = YTArray([3, 3, 3])

    yield operate_and_compare, a1, a2, operator.sub, answer1
    yield operate_and_compare, a2, a1, operator.sub, answer2
    yield operate_and_compare, a1, a3, operator.sub, answer1
    yield operate_and_compare, a3, a1, operator.sub, answer2
    yield operate_and_compare, a1, a2, np.subtract, answer1
    yield operate_and_compare, a2, a1, np.subtract, answer2
    yield operate_and_compare, a1, a3, np.subtract, answer1
    yield operate_and_compare, a3, a1, np.subtract, answer2

    # Catch the different dimensions error
    a1 = YTArray([1, 2, 3], 'm')
    a2 = YTArray([4, 5, 6], 'kg')

    yield assert_raises, YTUnitOperationError, operator.sub, a1, a2
    yield assert_raises, YTUnitOperationError, operator.isub, a1, a2

    # subtracting with zero is allowed irrespective of the units
    zeros = np.zeros(3)
    zeros_yta_dimless = YTArray(zeros, 'dimensionless')
    zeros_yta_length = YTArray(zeros, 'm')
    zeros_yta_mass = YTArray(zeros, 'kg')
    operands = [0, YTQuantity(0), YTQuantity(0, 'kg'), zeros, zeros_yta_dimless,
                zeros_yta_length, zeros_yta_mass]

    for op in [operator.sub, np.subtract]:
        for operand in operands:
            yield operate_and_compare, a1, operand, op, a1
            yield operate_and_compare, operand, a1, op, -a1
            yield operate_and_compare, 4*m, operand, op, 4*m
            yield operate_and_compare, operand, 4*m, op, -4*m

def test_multiplication():
    """
    Test multiplication of two YTArrays

    """

    # Same units
    a1 = YTArray([1, 2, 3], 'cm')
    a2 = YTArray([4, 5, 6], 'cm')
    a3 = [4*cm, 5*cm, 6*cm]
    answer = YTArray([4, 10, 18], 'cm**2')

    yield operate_and_compare, a1, a2, operator.mul, answer
    yield operate_and_compare, a2, a1, operator.mul, answer
    yield operate_and_compare, a1, a3, operator.mul, answer
    yield operate_and_compare, a3, a1, operator.mul, answer
    yield operate_and_compare, a1, a2, np.multiply, answer
    yield operate_and_compare, a2, a1, np.multiply, answer
    yield operate_and_compare, a1, a3, np.multiply, answer
    yield operate_and_compare, a3, a1, np.multiply, answer

    # different units, same dimension
    a1 = YTArray([1, 2, 3], 'cm')
    a2 = YTArray([4, 5, 6], 'm')
    a3 = [4*m, 5*m, 6*m]
    answer1 = YTArray([400, 1000, 1800], 'cm**2')
    answer2 = YTArray([.04, .10, .18], 'm**2')
    answer3 = YTArray([4, 10, 18], 'cm*m')

    yield operate_and_compare, a1, a2, operator.mul, answer1
    yield operate_and_compare, a2, a1, operator.mul, answer2
    yield operate_and_compare, a1, a3, operator.mul, answer1
    yield operate_and_compare, a3, a1, operator.mul, answer2
    yield operate_and_compare, a1, a2, np.multiply, answer3
    yield operate_and_compare, a2, a1, np.multiply, answer3
    yield operate_and_compare, a1, a3, np.multiply, answer3
    yield operate_and_compare, a3, a1, np.multiply, answer3

    # different dimensions
    a1 = YTArray([1, 2, 3], 'cm')
    a2 = YTArray([4, 5, 6], 'g')
    a3 = [4*g, 5*g, 6*g]
    answer = YTArray([4, 10, 18], 'cm*g')

    yield operate_and_compare, a1, a2, operator.mul, answer
    yield operate_and_compare, a2, a1, operator.mul, answer
    yield operate_and_compare, a1, a3, operator.mul, answer
    yield operate_and_compare, a3, a1, operator.mul, answer
    yield operate_and_compare, a1, a2, np.multiply, answer
    yield operate_and_compare, a2, a1, np.multiply, answer
    yield operate_and_compare, a1, a3, np.multiply, answer
    yield operate_and_compare, a3, a1, np.multiply, answer

    # One dimensionless, one unitful
    a1 = YTArray([1, 2, 3], 'cm')
    a2 = array([4, 5, 6])
    a3 = [4, 5, 6]
    answer = YTArray([4, 10, 18], 'cm')

    yield operate_and_compare, a1, a2, operator.mul, answer
    yield operate_and_compare, a2, a1, operator.mul, answer
    yield operate_and_compare, a1, a3, operator.mul, answer
    yield operate_and_compare, a3, a1, operator.mul, answer
    yield operate_and_compare, a1, a2, np.multiply, answer
    yield operate_and_compare, a2, a1, np.multiply, answer
    yield operate_and_compare, a1, a3, np.multiply, answer
    yield operate_and_compare, a3, a1, np.multiply, answer

    # Both dimensionless quantities
    a1 = YTArray([1, 2, 3])
    a2 = array([4, 5, 6])
    a3 = [4, 5, 6]
    answer = YTArray([4, 10, 18])

    yield operate_and_compare, a1, a2, operator.mul, answer
    yield operate_and_compare, a2, a1, operator.mul, answer
    yield operate_and_compare, a1, a3, operator.mul, answer
    yield operate_and_compare, a3, a1, operator.mul, answer
    yield operate_and_compare, a1, a2, np.multiply, answer
    yield operate_and_compare, a2, a1, np.multiply, answer
    yield operate_and_compare, a1, a3, np.multiply, answer
    yield operate_and_compare, a3, a1, np.multiply, answer


def test_division():
    """
    Test multiplication of two YTArrays

    """

    # Same units
    a1 = YTArray([1., 2., 3.], 'cm')
    a2 = YTArray([4., 5., 6.], 'cm')
    a3 = [4*cm, 5*cm, 6*cm]
    answer1 = YTArray([0.25, 0.4, 0.5])
    answer2 = YTArray([4, 2.5, 2])
    if "div" in dir(operator):
        op = operator.div
    else:
        op = operator.truediv

    yield operate_and_compare, a1, a2, op, answer1
    yield operate_and_compare, a2, a1, op, answer2
    yield operate_and_compare, a1, a3, op, answer1
    yield operate_and_compare, a3, a1, op, answer2
    yield operate_and_compare, a1, a2, np.divide, answer1
    yield operate_and_compare, a2, a1, np.divide, answer2
    yield operate_and_compare, a1, a3, np.divide, answer1
    yield operate_and_compare, a3, a1, np.divide, answer2

    # different units, same dimension
    a1 = YTArray([1., 2., 3.], 'cm')
    a2 = YTArray([4., 5., 6.], 'm')
    a3 = [4*m, 5*m, 6*m]
    answer1 = YTArray([.0025, .004, .005])
    answer2 = YTArray([400, 250, 200])
    answer3 = YTArray([0.25, 0.4, 0.5], 'cm/m')
    answer4 = YTArray([4.0, 2.5, 2.0], 'm/cm')

    yield operate_and_compare, a1, a2, op, answer1
    yield operate_and_compare, a2, a1, op, answer2
    yield operate_and_compare, a1, a3, op, answer1
    yield operate_and_compare, a3, a1, op, answer2
    yield operate_and_compare, a1, a2, np.divide, answer3
    yield operate_and_compare, a2, a1, np.divide, answer4
    yield operate_and_compare, a1, a3, np.divide, answer3
    yield operate_and_compare, a3, a1, np.divide, answer4

    # different dimensions
    a1 = YTArray([1., 2., 3.], 'cm')
    a2 = YTArray([4., 5., 6.], 'g')
    a3 = [4*g, 5*g, 6*g]
    answer1 = YTArray([0.25, 0.4, 0.5], 'cm/g')
    answer2 = YTArray([4, 2.5, 2], 'g/cm')

    yield operate_and_compare, a1, a2, op, answer1
    yield operate_and_compare, a2, a1, op, answer2
    yield operate_and_compare, a1, a3, op, answer1
    yield operate_and_compare, a3, a1, op, answer2
    yield operate_and_compare, a1, a2, np.divide, answer1
    yield operate_and_compare, a2, a1, np.divide, answer2
    yield operate_and_compare, a1, a3, np.divide, answer1
    yield operate_and_compare, a3, a1, np.divide, answer2

    # One dimensionless, one unitful
    a1 = YTArray([1., 2., 3.], 'cm')
    a2 = array([4., 5., 6.])
    a3 = [4, 5, 6]
    answer1 = YTArray([0.25, 0.4, 0.5], 'cm')
    answer2 = YTArray([4, 2.5, 2], '1/cm')

    yield operate_and_compare, a1, a2, op, answer1
    yield operate_and_compare, a2, a1, op, answer2
    yield operate_and_compare, a1, a3, op, answer1
    yield operate_and_compare, a3, a1, op, answer2
    yield operate_and_compare, a1, a2, np.divide, answer1
    yield operate_and_compare, a2, a1, np.divide, answer2
    yield operate_and_compare, a1, a3, np.divide, answer1
    yield operate_and_compare, a3, a1, np.divide, answer2

    # Both dimensionless quantities
    a1 = YTArray([1., 2., 3.])
    a2 = array([4., 5., 6.])
    a3 = [4, 5, 6]
    answer1 = YTArray([0.25, 0.4, 0.5])
    answer2 = YTArray([4, 2.5, 2])

    yield operate_and_compare, a1, a2, op, answer1
    yield operate_and_compare, a2, a1, op, answer2
    yield operate_and_compare, a1, a3, op, answer1
    yield operate_and_compare, a3, a1, op, answer2
    yield operate_and_compare, a1, a3, np.divide, answer1
    yield operate_and_compare, a3, a1, np.divide, answer2
    yield operate_and_compare, a1, a3, np.divide, answer1
    yield operate_and_compare, a3, a1, np.divide, answer2


def test_power():
    """
    Test power operator ensure units are correct.

    """

    from yt.units import cm

    cm_arr = np.array([1.0, 1.0]) * cm

    assert_equal, cm**3, YTQuantity(1, 'cm**3')
    assert_equal, np.power(cm, 3), YTQuantity(1, 'cm**3')
    assert_equal, cm**YTQuantity(3), YTQuantity(1, 'cm**3')
    assert_raises, YTUnitOperationError, np.power, cm, YTQuantity(3, 'g')

    assert_equal, cm_arr**3, YTArray([1, 1], 'cm**3')
    assert_equal, np.power(cm_arr, 3), YTArray([1, 1], 'cm**3')
    assert_equal, cm_arr**YTQuantity(3), YTArray([1, 1], 'cm**3')
    assert_raises, YTUnitOperationError, np.power, cm_arr, YTQuantity(3, 'g')


def test_comparisons():
    """
    Test numpy ufunc comparison operators for unit consistency.

    """
    from yt.units.yt_array import YTArray

    a1 = YTArray([1, 2, 3], 'cm')
    a2 = YTArray([2, 1, 3], 'cm')
    a3 = YTArray([.02, .01, .03], 'm')
    dimless = np.array([2,1,3])

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
    for op, answer in zip(ops, answers):
        yield operate_and_compare, a1, dimless, op, answer

    for op in ops:
        yield assert_raises, YTUfuncUnitError, op, a1, a3

    for op, answer in zip(ops, answers):
        yield operate_and_compare, a1, a3.in_units('cm'), op, answer
    
    # Check that comparisons with dimensionless quantities work in both directions.
    yield operate_and_compare, a3, dimless, np.less, [True, True, True]
    yield operate_and_compare, dimless, a3, np.less, [False, False, False]
    yield assert_equal, a1 < 2, [True, False, False]
    yield assert_equal, a1 < 2, np.less(a1, 2)
    yield assert_equal, 2 < a1, [False, False, True]
    yield assert_equal, 2 < a1, np.less(2, a1)


def test_unit_conversions():
    """
    Test operations that convert to different units or cast to ndarray

    """
    from yt.units.yt_array import YTQuantity
    from yt.units.unit_object import Unit

    km = YTQuantity(1, 'km')
    km_in_cm = km.in_units('cm')
    cm_unit = Unit('cm')
    kpc_unit = Unit('kpc')

    yield assert_equal, km_in_cm, km
    yield assert_equal, km_in_cm.in_cgs(), 1e5
    yield assert_equal, km_in_cm.in_mks(), 1e3
    yield assert_equal, km_in_cm.units, cm_unit

    km_view = km.ndarray_view()
    km.convert_to_units('cm')
    assert_true(km_view.base is km.base)

    yield assert_equal, km, YTQuantity(1, 'km')
    yield assert_equal, km.in_cgs(), 1e5
    yield assert_equal, km.in_mks(), 1e3
    yield assert_equal, km.units, cm_unit

    km.convert_to_units('kpc')
    assert_true(km_view.base is km.base)

    yield assert_array_almost_equal_nulp, km, YTQuantity(1, 'km')
    yield assert_array_almost_equal_nulp, km.in_cgs(), YTQuantity(1e5, 'cm')
    yield assert_array_almost_equal_nulp, km.in_mks(), YTQuantity(1e3, 'm')
    yield assert_equal, km.units, kpc_unit

    yield assert_isinstance, km.to_ndarray(), np.ndarray
    yield assert_isinstance, km.ndarray_view(), np.ndarray

    dyne = YTQuantity(1.0, 'dyne')

    yield assert_equal, dyne.in_cgs(), dyne
    yield assert_equal, dyne.in_cgs(), 1.0
    yield assert_equal, dyne.in_mks(), dyne
    yield assert_equal, dyne.in_mks(), 1e-5
    yield assert_equal, str(dyne.in_mks().units), 'kg*m/s**2'
    yield assert_equal, str(dyne.in_cgs().units), 'cm*g/s**2'

    em3 = YTQuantity(1.0, 'erg/m**3')

    yield assert_equal, em3.in_cgs(), em3
    yield assert_equal, em3.in_cgs(), 1e-6
    yield assert_equal, em3.in_mks(), em3
    yield assert_equal, em3.in_mks(), 1e-7
    yield assert_equal, str(em3.in_mks().units), 'kg/(m*s**2)'
    yield assert_equal, str(em3.in_cgs().units), 'g/(cm*s**2)'

    em3_converted = YTQuantity(1545436840.386756, 'Msun/(Myr**2*kpc)')
    yield assert_equal, em3.in_base(unit_system="galactic"), em3
    yield assert_array_almost_equal, em3.in_base(unit_system="galactic"), em3_converted
    yield assert_equal, str(em3.in_base(unit_system="galactic").units), 'Msun/(Myr**2*kpc)'

    dimless = YTQuantity(1.0, "")
    yield assert_equal, dimless.in_cgs(), dimless
    yield assert_equal, dimless.in_cgs(), 1.0
    yield assert_equal, dimless.in_mks(), dimless
    yield assert_equal, dimless.in_mks(), 1.0
    yield assert_equal, str(dimless.in_cgs().units), "dimensionless"

def test_temperature_conversions():
    """
    Test conversions between various supported temperatue scales.

    Also ensure we only allow compound units with temperature
    scales that have a proper zero point.

    """
    from yt.units.unit_object import InvalidUnitOperation

    km = YTQuantity(1, 'km')
    balmy = YTQuantity(300, 'K')
    balmy_F = YTQuantity(80.33, 'degF')
    balmy_C = YTQuantity(26.85, 'degC')
    balmy_R = YTQuantity(540, 'R')

    assert_array_almost_equal(balmy.in_units('degF'), balmy_F)
    assert_array_almost_equal(balmy.in_units('degC'), balmy_C)
    assert_array_almost_equal(balmy.in_units('R'), balmy_R)

    balmy_view = balmy.ndarray_view()

    balmy.convert_to_units('degF')
    yield assert_true, balmy_view.base is balmy.base
    yield assert_array_almost_equal, np.array(balmy), np.array(balmy_F)

    balmy.convert_to_units('degC')
    yield assert_true, balmy_view.base is balmy.base
    yield assert_array_almost_equal, np.array(balmy), np.array(balmy_C)

    balmy.convert_to_units('R')
    yield assert_true, balmy_view.base is balmy.base
    yield assert_array_almost_equal, np.array(balmy), np.array(balmy_R)

    balmy.convert_to_units('degF')
    yield assert_true, balmy_view.base is balmy.base
    yield assert_array_almost_equal, np.array(balmy), np.array(balmy_F)

    yield assert_raises, InvalidUnitOperation, np.multiply, balmy, km

    # Does CGS conversion from F to K work?
    yield assert_array_almost_equal, balmy.in_cgs(), YTQuantity(300, 'K')


def test_yt_array_yt_quantity_ops():
    """
    Test operations that combine YTArray and YTQuantity
    """
    a = YTArray(range(10), 'cm')
    b = YTQuantity(5, 'g')

    assert_isinstance(a*b, YTArray)
    assert_isinstance(b*a, YTArray)

    assert_isinstance(a/b, YTArray)
    assert_isinstance(b/a, YTArray)

    assert_isinstance(a*a, YTArray)
    assert_isinstance(a/a, YTArray)

    assert_isinstance(b*b, YTQuantity)
    assert_isinstance(b/b, YTQuantity)


def test_selecting():
    """
    Test slicing of two YTArrays

    """
    a = YTArray(range(10), 'cm')
    a_slice = a[:3]
    a_fancy_index = a[[1, 1, 3, 5]]
    a_array_fancy_index = a[array([[1, 1], [3, 5]])]
    a_boolean_index = a[a > 5]
    a_selection = a[0]

    yield assert_array_equal, a_slice, YTArray([0, 1, 2], 'cm')
    yield assert_array_equal, a_fancy_index, YTArray([1, 1, 3, 5], 'cm')
    yield assert_array_equal, a_array_fancy_index, \
        YTArray([[1, 1, ], [3, 5]], 'cm')
    yield assert_array_equal, a_boolean_index, YTArray([6, 7, 8, 9], 'cm')
    yield assert_isinstance, a_selection, YTQuantity

    # .base points to the original array for a numpy view.  If it is not a
    # view, .base is None.
    yield assert_true, a_slice.base is a


def test_iteration():
    """
    Test that iterating over a YTArray returns a sequence of YTQuantity insances
    """
    a = np.arange(3)
    b = YTArray(np.arange(3), 'cm')
    for ia, ib, in zip(a, b):
        yield assert_equal, ia, ib.value
        yield assert_equal, ib.units, b.units


def test_fix_length():
    """
    Test fixing the length of an array. Used in spheres and other data objects
    """
    ds = fake_random_ds(64, nprocs=1, length_unit=10)
    length = ds.quan(1.0, 'code_length')
    new_length = fix_length(length, ds=ds)
    yield assert_equal, YTQuantity(10, 'cm'), new_length

def test_code_unit_combinations():
    """
    Test comparing code units coming from different datasets
    """
    ds1 = fake_random_ds(64, nprocs=1, length_unit=1)
    ds2 = fake_random_ds(64, nprocs=1, length_unit=10)

    q1 = ds1.quan(1, 'code_length')
    q2 = ds2.quan(1, 'code_length')

    assert_equal(10*q1, q2)
    assert_equal(q1/q2, 0.1)
    assert_true(q1 < q2)
    assert_true(q2 > q1)
    assert_true(not bool(q1 > q2))
    assert_true(not bool(q2 < q1))
    assert_true(q1 != q2)
    assert_true(not bool(q1 == q2))

    assert_equal((q1 + q2).in_cgs().value, 11)
    assert_equal((q2 + q1).in_cgs().value, 11)
    assert_equal((q1 - q2).in_cgs().value, -9)
    assert_equal((q2 - q1).in_cgs().value, 9)

def test_ytarray_pickle():
    ds = fake_random_ds(64, nprocs=1)
    test_data = [ds.quan(12.0, 'code_length'),
                 ds.arr([1, 2, 3], 'code_length')]

    for data in test_data:
        tempf = tempfile.NamedTemporaryFile(delete=False)
        pickle.dump(data, tempf)
        tempf.close()

        with open(tempf.name, "rb") as fname:
            loaded_data = pickle.load(fname)
        os.unlink(tempf.name)

        yield assert_array_equal, data, loaded_data
        yield assert_equal, data.units, loaded_data.units
        yield assert_array_equal, array(data.in_cgs()), \
            array(loaded_data.in_cgs())
        yield assert_equal, float(data.units.base_value), \
            float(loaded_data.units.base_value)


def test_copy():
    quan = YTQuantity(1, 'g')
    arr = YTArray([1, 2, 3], 'cm')

    yield assert_equal, copy.copy(quan), quan
    yield assert_array_equal, copy.copy(arr), arr

    yield assert_equal,  copy.deepcopy(quan), quan
    yield assert_array_equal, copy.deepcopy(arr), arr

    yield assert_equal, quan.copy(), quan
    yield assert_array_equal, arr.copy(), arr

    yield assert_equal, np.copy(quan), quan
    yield assert_array_equal, np.copy(arr), arr


def unary_ufunc_comparison(ufunc, a):
    out = a.copy()
    a_array = a.to_ndarray()
    if ufunc in (np.isreal, np.iscomplex, ):
        # According to the numpy docs, these two explicitly do not do
        # in-place copies.
        ret = ufunc(a)
        assert_true(not hasattr(ret, 'units'))
        assert_array_equal(ret, ufunc(a))
    elif ufunc in (np.exp, np.exp2, np.log, np.log2, np.log10, np.expm1,
                   np.log1p, np.sin, np.cos, np.tan, np.arcsin, np.arccos,
                   np.arctan, np.sinh, np.cosh, np.tanh, np.arccosh,
                   np.arcsinh, np.arctanh, np.deg2rad, np.rad2deg,
                   np.isfinite, np.isinf, np.isnan, np.signbit, np.sign,
                   np.rint, np.logical_not):
        # These operations should return identical results compared to numpy.

        try:
            ret = ufunc(a, out=out)
        except YTUnitOperationError:
            assert_true(ufunc in (np.deg2rad, np.rad2deg))
            ret = ufunc(YTArray(a, '1'))

        assert_array_equal(ret, out)
        assert_array_equal(ret, ufunc(a_array))
        # In-place copies do not drop units.
        assert_true(hasattr(out, 'units'))
        assert_true(not hasattr(ret, 'units'))
    elif ufunc in (np.absolute, np.fabs, np.conjugate, np.floor, np.ceil,
                   np.trunc, np.negative, np.spacing):
        ret = ufunc(a, out=out)

        assert_array_equal(ret, out)
        assert_array_equal(ret.to_ndarray(), ufunc(a_array))
        assert_true(ret.units == out.units)
    elif ufunc in (np.ones_like, np.square, np.sqrt, np.reciprocal):
        if ufunc is np.ones_like:
            ret = ufunc(a)
        else:
            ret = ufunc(a, out=out)
            assert_array_equal(ret, out)

        assert_array_equal(ret.to_ndarray(), ufunc(a_array))
        if ufunc is np.square:
            assert_true(out.units == a.units**2)
            assert_true(ret.units == a.units**2)
        elif ufunc is np.sqrt:
            assert_true(out.units == a.units**0.5)
            assert_true(ret.units == a.units**0.5)
        elif ufunc is np.reciprocal:
            assert_true(out.units == a.units**-1)
            assert_true(ret.units == a.units**-1)
    elif ufunc is np.modf:
        ret1, ret2 = ufunc(a)
        npret1, npret2 = ufunc(a_array)

        assert_array_equal(ret1.to_ndarray(), npret1)
        assert_array_equal(ret2.to_ndarray(), npret2)
    elif ufunc is np.frexp:
        ret1, ret2 = ufunc(a)
        npret1, npret2 = ufunc(a_array)

        assert_array_equal(ret1, npret1)
        assert_array_equal(ret2, npret2)
    elif ufunc is np.invert:
        assert_raises(TypeError, ufunc, a)
    else:
        # There shouldn't be any untested ufuncs.
        assert_true(False)


def binary_ufunc_comparison(ufunc, a, b):
    out = a.copy()
    if ufunc in (np.add, np.subtract, np.remainder, np.fmod, np.mod,
                 np.arctan2, np.hypot, np.greater, np.greater_equal, np.less,
                 np.less_equal, np.equal, np.not_equal, np.logical_and,
                 np.logical_or, np.logical_xor, np.maximum, np.minimum,
                 np.fmax, np.fmin, np.nextafter):
        if a.units != b.units and a.units.dimensions == b.units.dimensions:
            assert_raises(YTUfuncUnitError, ufunc, a, b)
            return
        elif a.units != b.units:
            assert_raises(YTUnitOperationError, ufunc, a, b)
            return
    if ufunc in (np.bitwise_and, np.bitwise_or, np.bitwise_xor,
                 np.left_shift, np.right_shift, np.ldexp):
        assert_raises(TypeError, ufunc, a, b)
        return

    ret = ufunc(a, b, out=out)

    if ufunc is np.multiply:
        assert_true(ret.units == a.units*b.units)
    elif ufunc in (np.divide, np.true_divide, np.arctan2):
        assert_true(ret.units.dimensions == (a.units/b.units).dimensions)
    elif ufunc in (np.greater, np.greater_equal, np.less, np.less_equal,
                   np.not_equal, np.equal, np.logical_and, np.logical_or,
                   np.logical_xor):
        assert_true(not isinstance(ret, YTArray) and
                    isinstance(ret, np.ndarray))
    assert_array_equal(ret, out)
    if (ufunc in (np.divide, np.true_divide, np.arctan2) and
        (a.units.dimensions == b.units.dimensions)):
        assert_array_almost_equal(
            np.array(ret), ufunc(np.array(a.in_cgs()), np.array(b.in_cgs())))
    else:
        assert_array_almost_equal(np.array(ret), ufunc(np.array(a), np.array(b)))


def test_ufuncs():
    for ufunc in unary_operators:
        yield unary_ufunc_comparison, ufunc, YTArray([.3, .4, .5], 'cm')
        yield unary_ufunc_comparison, ufunc, YTArray([12, 23, 47], 'g')
        yield unary_ufunc_comparison, ufunc, YTArray([2, 4, -6], 'erg/m**3')

    for ufunc in binary_operators:

        # arr**arr is undefined for arrays with units because
        # each element of the result would have different units.
        if ufunc is np.power:
            a = YTArray([.3, .4, .5], 'cm')
            b = YTArray([.1, .2, .3], 'dimensionless')
            c = np.array(b)
            d = YTArray([1., 2., 3.], 'g')
            yield binary_ufunc_comparison, ufunc, a, b
            yield binary_ufunc_comparison, ufunc, a, c
            assert_raises(YTUnitOperationError, ufunc, a, d)
            continue

        a = YTArray([.3, .4, .5], 'cm')
        b = YTArray([.1, .2, .3], 'cm')
        c = YTArray([.1, .2, .3], 'm')
        d = YTArray([.1, .2, .3], 'g')
        e = YTArray([.1, .2, .3], 'erg/m**3')

        for pair in itertools.product([a, b, c, d, e], repeat=2):
            yield binary_ufunc_comparison, ufunc, pair[0], pair[1]


def test_convenience():

    arr = YTArray([1, 2, 3], 'cm')

    yield assert_equal, arr.unit_quantity, YTQuantity(1, 'cm')
    yield assert_equal, arr.uq, YTQuantity(1, 'cm')
    yield assert_isinstance, arr.unit_quantity, YTQuantity
    yield assert_isinstance, arr.uq, YTQuantity

    yield assert_array_equal, arr.unit_array, YTArray(np.ones_like(arr), 'cm')
    yield assert_array_equal, arr.ua, YTArray(np.ones_like(arr), 'cm')
    yield assert_isinstance, arr.unit_array, YTArray
    yield assert_isinstance, arr.ua, YTArray

    yield assert_array_equal, arr.ndview, arr.view(np.ndarray)
    yield assert_array_equal, arr.d, arr.view(np.ndarray)
    yield assert_true, arr.ndview.base is arr.base
    yield assert_true, arr.d.base is arr.base

    yield assert_array_equal, arr.value, np.array(arr)
    yield assert_array_equal, arr.v, np.array(arr)


def test_registry_association():
    ds = fake_random_ds(64, nprocs=1, length_unit=10)
    a = ds.quan(3, 'cm')
    b = YTQuantity(4, 'm')
    c = ds.quan(6, '')
    d = 5

    yield assert_equal, id(a.units.registry), id(ds.unit_registry)

    def binary_op_registry_comparison(op):
        e = op(a, b)
        f = op(b, a)
        g = op(c, d)
        h = op(d, c)

        assert_equal(id(e.units.registry), id(ds.unit_registry))
        assert_equal(id(f.units.registry), id(b.units.registry))
        assert_equal(id(g.units.registry), id(h.units.registry))
        assert_equal(id(g.units.registry), id(ds.unit_registry))

    def unary_op_registry_comparison(op):
        c = op(a)
        d = op(b)

        assert_equal(id(c.units.registry), id(ds.unit_registry))
        assert_equal(id(d.units.registry), id(b.units.registry))

    binary_ops = [operator.add, operator.sub, operator.mul, 
                  operator.truediv]
    if hasattr(operator, "div"):
        binary_ops.append(operator.div)
    for op in binary_ops:
        yield binary_op_registry_comparison, op

    for op in [operator.abs, operator.neg, operator.pos]:
        yield unary_op_registry_comparison, op

@requires_module("astropy")
def test_astropy():
    from yt.utilities.on_demand_imports import _astropy

    ap_arr = np.arange(10)*_astropy.units.km/_astropy.units.hr
    yt_arr = YTArray(np.arange(10), "km/hr")
    yt_arr2 = YTArray.from_astropy(ap_arr)

    ap_quan = 10.*_astropy.units.Msun**0.5/(_astropy.units.kpc**3)
    yt_quan = YTQuantity(10., "sqrt(Msun)/kpc**3")
    yt_quan2 = YTQuantity.from_astropy(ap_quan)

    yield assert_array_equal, ap_arr, yt_arr.to_astropy()
    yield assert_array_equal, yt_arr, YTArray.from_astropy(ap_arr)
    yield assert_array_equal, yt_arr, yt_arr2

    yield assert_equal, ap_quan, yt_quan.to_astropy()
    yield assert_equal, yt_quan, YTQuantity.from_astropy(ap_quan)
    yield assert_equal, yt_quan, yt_quan2

    yield assert_array_equal, yt_arr, YTArray.from_astropy(yt_arr.to_astropy())
    yield assert_equal, yt_quan, YTQuantity.from_astropy(yt_quan.to_astropy())

@requires_module("pint")
def test_pint():
    from pint import UnitRegistry

    ureg = UnitRegistry()
    
    p_arr = np.arange(10)*ureg.km/ureg.hr
    yt_arr = YTArray(np.arange(10), "km/hr")
    yt_arr2 = YTArray.from_pint(p_arr)

    p_quan = 10.*ureg.g**0.5/(ureg.mm**3)
    yt_quan = YTQuantity(10., "sqrt(g)/mm**3")
    yt_quan2 = YTQuantity.from_pint(p_quan)

    yield assert_array_equal, p_arr, yt_arr.to_pint()
    assert_equal(p_quan, yt_quan.to_pint())
    yield assert_array_equal, yt_arr, YTArray.from_pint(p_arr)
    yield assert_array_equal, yt_arr, yt_arr2

    yield assert_equal, p_quan.magnitude, yt_quan.to_pint().magnitude
    assert_equal(p_quan, yt_quan.to_pint())
    yield assert_equal, yt_quan, YTQuantity.from_pint(p_quan)
    yield assert_equal, yt_quan, yt_quan2

    yield assert_array_equal, yt_arr, YTArray.from_pint(yt_arr.to_pint())
    yield assert_equal, yt_quan, YTQuantity.from_pint(yt_quan.to_pint())

def test_subclass():

    class YTASubclass(YTArray):
        pass

    a = YTASubclass([4, 5, 6], 'g')
    b = YTASubclass([7, 8, 9], 'kg')
    nu = YTASubclass([10, 11, 12], '')
    nda = np.array([3, 4, 5])
    yta = YTArray([6, 7, 8], 'mg')
    loq = [YTQuantity(6, 'mg'), YTQuantity(7, 'mg'), YTQuantity(8, 'mg')]
    ytq = YTQuantity(4, 'cm')
    ndf = np.float64(3)

    def op_comparison(op, inst1, inst2, compare_class):
        assert_isinstance(op(inst1, inst2), compare_class)
        assert_isinstance(op(inst2, inst1), compare_class)

    ops = [operator.mul, operator.truediv]
    if hasattr(operator, "div"):
        ops.append(operator.div)
    for op in ops:
        for inst in (b, ytq, ndf, yta, nda, loq):
            op_comparison(op, a, inst, YTASubclass)

        op_comparison(op, ytq, nda, YTArray)
        op_comparison(op, ytq, yta, YTArray)

    for op in (operator.add, operator.sub):
        op_comparison(op, nu, nda, YTASubclass)
        op_comparison(op, a, b, YTASubclass)
        op_comparison(op, a, yta, YTASubclass)
        op_comparison(op, a, loq, YTASubclass)

    assert_isinstance(a[0], YTQuantity)
    assert_isinstance(a[:], YTASubclass)
    assert_isinstance(a[:2], YTASubclass)
    assert_isinstance(YTASubclass(yta), YTASubclass)
    
def test_h5_io():
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    ds = fake_random_ds(64, nprocs=1, length_unit=10)

    warr = ds.arr(np.random.random((256, 256)), 'code_length')

    warr.write_hdf5('test.h5')

    iarr = YTArray.from_hdf5('test.h5')

    yield assert_equal, warr, iarr
    yield assert_equal, warr.units.registry['code_length'], iarr.units.registry['code_length']

    warr.write_hdf5('test.h5', dataset_name="test_dset", group_name='/arrays/test_group')

    giarr = YTArray.from_hdf5('test.h5', dataset_name="test_dset", group_name='/arrays/test_group')

    yield assert_equal, warr, giarr

    os.chdir(curdir)
    shutil.rmtree(tmpdir)

def test_equivalencies():
    from yt.utilities.physical_constants import clight, mp, kboltz, hcgs, mh, me, \
        mass_sun_cgs, G, stefan_boltzmann_constant_cgs
    import yt.units as u

    # Mass-energy

    E = mp.to_equivalent("keV","mass_energy")
    yield assert_equal, E, mp*clight*clight
    yield assert_allclose_units, mp, E.to_equivalent("g", "mass_energy")

    # Thermal

    T = YTQuantity(1.0e8,"K")
    E = T.to_equivalent("W*hr","thermal")
    yield assert_equal, E, (kboltz*T).in_units("W*hr")
    yield assert_allclose_units, T, E.to_equivalent("K", "thermal")

    # Spectral

    l = YTQuantity(4000.,"angstrom")
    nu = l.to_equivalent("Hz","spectral")
    yield assert_equal, nu, clight/l
    E = hcgs*nu
    l2 = E.to_equivalent("angstrom", "spectral")
    yield assert_allclose_units, l, l2
    nu2 = clight/l2.in_units("cm")
    yield assert_allclose_units, nu, nu2
    E2 = nu2.to_equivalent("keV", "spectral")
    yield assert_allclose_units, E2, E.in_units("keV")

    # Sound-speed

    mu = 0.6
    gg = 5./3.
    c_s = T.to_equivalent("km/s","sound_speed")
    yield assert_equal, c_s, np.sqrt(gg*kboltz*T/(mu*mh))
    yield assert_allclose_units, T, c_s.to_equivalent("K","sound_speed")

    mu = 0.5
    gg = 4./3.
    c_s = T.to_equivalent("km/s","sound_speed", mu=mu, gamma=gg)
    yield assert_equal, c_s, np.sqrt(gg*kboltz*T/(mu*mh))
    yield assert_allclose_units, T, c_s.to_equivalent("K","sound_speed",
                                                    mu=mu, gamma=gg)

    # Lorentz

    v = 0.8*clight
    g = v.to_equivalent("dimensionless","lorentz")
    g2 = YTQuantity(1./np.sqrt(1.-0.8*0.8), "dimensionless")
    yield assert_allclose_units, g, g2
    v2 = g2.to_equivalent("mile/hr", "lorentz")
    yield assert_allclose_units, v2, v.in_units("mile/hr")

    # Schwarzschild

    R = mass_sun_cgs.to_equivalent("kpc","schwarzschild")
    yield assert_equal, R.in_cgs(), 2*G*mass_sun_cgs/(clight*clight)
    yield assert_allclose_units, mass_sun_cgs, R.to_equivalent("g", "schwarzschild")

    # Compton

    l = me.to_equivalent("angstrom","compton")
    yield assert_equal, l, hcgs/(me*clight)
    yield assert_allclose_units, me, l.to_equivalent("g", "compton")

    # Number density

    rho = mp/u.cm**3

    n = rho.to_equivalent("cm**-3","number_density")
    yield assert_equal, n, rho/(mh*0.6)
    yield assert_allclose_units, rho, n.to_equivalent("g/cm**3","number_density")

    n = rho.to_equivalent("cm**-3","number_density", mu=0.75)
    yield assert_equal, n, rho/(mh*0.75)
    yield assert_allclose_units, rho, n.to_equivalent("g/cm**3","number_density", mu=0.75)

    # Effective temperature

    T = YTQuantity(1.0e4, "K")
    F = T.to_equivalent("erg/s/cm**2","effective_temperature")
    yield assert_equal, F, stefan_boltzmann_constant_cgs*T**4
    yield assert_allclose_units, T, F.to_equivalent("K", "effective_temperature")

def test_electromagnetic():
    from yt.units.dimensions import charge_mks, pressure, current_cgs, \
        magnetic_field_mks, magnetic_field_cgs, power
    from yt.utilities.physical_constants import mu_0, qp
    from yt.utilities.physical_ratios import speed_of_light_cm_per_s

    # Various tests of SI and CGS electromagnetic units

    qp_mks = qp.to_equivalent("C", "SI")
    yield assert_equal, qp_mks.units.dimensions, charge_mks
    yield assert_array_almost_equal, qp_mks.v, 10.0*qp.v/speed_of_light_cm_per_s

    qp_cgs = qp_mks.to_equivalent("esu", "CGS")
    yield assert_array_almost_equal, qp_cgs, qp
    yield assert_equal, qp_cgs.units.dimensions, qp.units.dimensions
    
    qp_mks_k = qp.to_equivalent("kC", "SI")
    yield assert_array_almost_equal, qp_mks_k.v, 1.0e-2*qp.v/speed_of_light_cm_per_s

    B = YTQuantity(1.0, "T")
    B_cgs = B.to_equivalent("gauss", "CGS")
    yield assert_equal, B.units.dimensions, magnetic_field_mks
    yield assert_equal, B_cgs.units.dimensions, magnetic_field_cgs
    yield assert_array_almost_equal, B_cgs, YTQuantity(1.0e4, "gauss")

    u_mks = B*B/(2*mu_0)
    yield assert_equal, u_mks.units.dimensions, pressure
    u_cgs = B_cgs*B_cgs/(8*np.pi)
    yield assert_equal, u_cgs.units.dimensions, pressure
    yield assert_array_almost_equal, u_mks.in_cgs(), u_cgs
    
    I = YTQuantity(1.0, "A")
    I_cgs = I.to_equivalent("statA", "CGS")
    yield assert_array_almost_equal, I_cgs, YTQuantity(0.1*speed_of_light_cm_per_s, "statA")
    yield assert_array_almost_equal, I_cgs.to_equivalent("mA", "SI"), I.in_units("mA")
    yield assert_equal, I_cgs.units.dimensions, current_cgs
    
    R = YTQuantity(1.0, "ohm")
    R_cgs = R.to_equivalent("statohm", "CGS")
    P_mks = I*I*R
    P_cgs = I_cgs*I_cgs*R_cgs
    yield assert_equal, P_mks.units.dimensions, power
    yield assert_equal, P_cgs.units.dimensions, power
    yield assert_array_almost_equal, P_cgs.in_cgs(), P_mks.in_cgs()
    yield assert_array_almost_equal, P_cgs.in_mks(), YTQuantity(1.0, "W")
    
    V = YTQuantity(1.0, "statV")
    V_mks = V.to_equivalent("V", "SI")
    yield assert_array_almost_equal, V_mks.v, 1.0e8*V.v/speed_of_light_cm_per_s

def test_ytarray_coercion():
    a = YTArray([1, 2, 3], 'cm')
    q = YTQuantity(3, 'cm')
    na = np.array([1, 2, 3])

    assert_isinstance(a*q, YTArray)
    assert_isinstance(q*na, YTArray)
    assert_isinstance(q*3, YTQuantity)
    assert_isinstance(q*np.float64(3), YTQuantity)
    assert_isinstance(q*np.array(3), YTQuantity)

def test_numpy_wrappers():
    a1 = YTArray([1, 2, 3], 'cm')
    a2 = YTArray([2, 3, 4, 5, 6], 'cm')
    catenate_answer = [1, 2, 3, 2, 3, 4, 5, 6]
    intersect_answer = [2, 3]
    union_answer = [1, 2, 3, 4, 5, 6]

    yield (assert_array_equal, YTArray(catenate_answer, 'cm'),
           uconcatenate((a1, a2)))
    yield assert_array_equal, catenate_answer, np.concatenate((a1, a2))

    yield (assert_array_equal, YTArray(intersect_answer, 'cm'),
           uintersect1d(a1, a2))
    yield assert_array_equal, intersect_answer, np.intersect1d(a1, a2)

    yield assert_array_equal, YTArray(union_answer, 'cm'), uunion1d(a1, a2)
    yield assert_array_equal, union_answer, np.union1d(a1, a2)

def test_dimensionless_conversion():
    a = YTQuantity(1, 'Zsun')
    b = a.in_units('Zsun')
    a.convert_to_units('Zsun')
    yield assert_true, a.units.base_value == metallicity_sun
    yield assert_true, b.units.base_value == metallicity_sun

def test_modified_unit_division():
    ds1 = fake_random_ds(64)
    ds2 = fake_random_ds(64)

    # this mocks comoving coordinates without going through the trouble
    # of setting up a fake cosmological dataset
    ds1.unit_registry.modify('m', 50)

    a = ds1.quan(3, 'm')
    b = ds2.quan(3, 'm')

    ret = a/b
    yield assert_true, ret == 0.5
    yield assert_true, ret.units.is_dimensionless
    yield assert_true, ret.units.base_value == 1.0

def test_load_and_save():
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    a = YTArray(np.random.random(10), "kpc")
    b = YTArray(np.random.random(10), "Msun")
    c = YTArray(np.random.random(10), "km/s")

    savetxt("arrays.dat", [a,b,c], delimiter=",")

    d, e = loadtxt("arrays.dat", usecols=(1,2), delimiter=",")

    yield assert_array_equal, b, d
    yield assert_array_equal, c, e

    os.chdir(curdir)
    shutil.rmtree(tmpdir)

def test_trig_ufunc_degrees():
    for ufunc in (np.sin, np.cos, np.tan):
        degree_values = np.random.random(10)*degree
        radian_values = degree_values.in_units('radian')
        assert_array_equal(ufunc(degree_values), ufunc(radian_values))

def test_builtin_sum():
    from yt.units import km

    arr = [1, 2, 3]*km
    assert_equal(sum(arr), 6*km)

def test_initialization_different_registries():
    from yt.testing import fake_random_ds

    ds1 = fake_random_ds(32, length_unit=1)
    ds2 = fake_random_ds(32, length_unit=3)

    l1 = ds1.quan(0.3, 'unitary')
    l2 = ds2.quan(l1, 'unitary')

    assert_almost_equal(float(l1.in_cgs()), 0.3)
    assert_almost_equal(float(l2.in_cgs()), 0.9)
    assert_almost_equal(float(ds1.quan(0.3, 'unitary').in_cgs()), 0.3)
    assert_almost_equal(float(ds2.quan(0.3, 'unitary').in_cgs()), 0.9)
