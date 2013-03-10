"""
YTArray class.

Author: Casey W. Stark <caseywstark@gmail.com>
Affiliation: UC Berkeley

Homepage: http://yt-project.org/
License:
  Copyright (C) 2013 Casey W. Stark.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""

import copy

import numpy as np

from numpy import add, subtract, multiply, divide, \
    negative, absolute, sqrt, square, power, reciprocal, ones_like, \
    isnan, isinf, cos, sin, log10, greater, equal, not_equal

from yt.utilities.units import Unit, UnitOperationError, dimensionless

def sqrt_unit(unit):
    return unit**0.5

def multiply_units(unit1, unit2):
    return unit1 * unit2

def preserve_unit(unit1, unit2):
    return unit1

def power_unit(unit, power):
    return unit**power

def square_unit(unit):
    return unit*unit

def divide_units(unit1, unit2):
    return unit1/unit2

def ones_like_units(unit):
    return unit

def isnan_unit(unit):
    return dimensionless

def isinf_unit(unit):
    return dimensionless

def negative_unit(unit):
    return unit

def absolute_unit(unit):
    return unit

def cos_unit(unit):
    return dimensionless

def sin_unit(unit):
    return dimensionless

def log10_unit(unit):
    return sympy.log(unit, 10)

def greater_unit(unit1, unit2):
    return dimensionless

def equal_unit(unit1, unit2):
    return dimensionless

def not_equal_unit(unit1, unit2):
    return dimensionless

class YTArray(np.ndarray):
    """

    """
    _ufunc_registry = {sqrt: sqrt_unit, multiply: multiply_units,
                       add: preserve_unit, subtract: preserve_unit,
                       power: power_unit, divide: divide_units,
                       square: square_unit, ones_like: ones_like_units,
                       isnan: isnan_unit, isinf: isinf_unit,
                       negative: negative_unit, absolute: absolute_unit,
                       cos: cos_unit, sin: sin_unit, log10: log10_unit,
                       greater: greater_unit, equal: equal_unit,
                       not_equal: not_equal_unit}

    def __new__(cls, input_array, input_units=None):
        if isinstance(input_array, YTArray):
            return input_array

        # Input array is an already formed ndarray instance
        # We first cast to be our class type
        obj = np.asarray(input_array).view(cls)

        # Check units type
        if input_units is None:
            # Nothing provided. Make dimensionless...
            units = Unit()
        elif isinstance(input_units, Unit):
            units = input_units
        else:
            # units kwarg set, but it's not a Unit object.
            # don't handle all the cases here, let the Unit class handle if
            # it's a str.
            units = Unit(input_units)

        # Attach the units
        obj.units = units

        return obj

    def __array_finalize__(self, obj):
        """

        """
        if obj is None:
            return
        self.units = getattr(obj, 'units', None)

    def __repr__(self):
        """

        """
        return super(YTArray, self).__repr__()+' '+self.units.__repr__()

    def __str__(self):
        """

        """
        return super(YTArray, self).__str__()+' '+self.units.__str__()

    #
    # Start unit conversion methods
    #

    def _unit_repr_check_same(self, units):
        """
        Takes a Unit object, or string of known unit symbol, and check that it
        is compatible with this quantity. Returns Unit object.

        """
        # let Unit() handle units arg if it's not already a Unit obj.
        if not isinstance(units, Unit):
            units = Unit(units)

        if not self.units.same_dimensions_as(units):
            raise UnitOperationError("Unit dimensionalities do not match. Tried to convert between %s (dim %s) and %s (dim %s)." % (self.units, self.units.dimensions, units, units.dimensions))

        return units

    def convert_to_units(self, units):
        """
        Convert the array and units to the given units.

        Parameters
        ----------
        units : Unit object or str
            The units you want to convert to.

        """
        new_units = self._unit_repr_check_same(units)
        conversion_factor = self.units.get_conversion_factor(new_units)

        self *= conversion_factor
        self.units = new_units

    def convert_to_cgs(self):
        """
        Convert the array and units to the equivalent cgs units.

        """
        return self.convert_to_units(self.units.get_cgs_equivalent())

    def in_units(self, units):
        """
        Creates a copy of this array with the data in the supplied units, and
        returns it.

        Parameters
        ----------
        units : Unit object or string
            The units you want to get a new quantity in.

        Returns
        -------
        YTArray

        """
        new_units = self._unit_repr_check_same(units)
        conversion_factor = self.units.get_conversion_factor(new_units)

        new_array = self * conversion_factor
        new_array.units = new_units

        return new_array

    def in_cgs(self):
        """
        Creates a copy of this array with the data in the equivalent cgs units,
        and returns it.

        Returns
        -------
        Quantity object with data converted to cgs and cgs units.

        """
        return self.in_units(self.units.get_cgs_equivalent())

    #
    # End unit conversion methods
    #

    #
    # Start operation methods
    #

    def __add__(self, right_object):
        """
        Add this ytarray to the object on the right of the `+` operator. Must
        check for the correct (same dimension) units.

        """
        # Make sure the other object is a YTArray before we use the `units`
        # attribute.
        if isinstance(right_object, YTArray):
            if not self.units.same_dimensions_as(right_object.units):
                raise UnitOperationError("You cannot add these YTArrays "
                    "because the unit dimensions do not match. `%s + %s` is "
                    "ill-defined." % (self.units, right_object.units))
        # If the other object is not a YTArray, the only valid case is adding
        # dimensionless things.
        else:
            if not self.units.is_dimensionless:
                raise UnitOperationError("You cannot add a pure number to a "
                    "dimensional YTArray. `%s + %s` is ill-defined."
                    % (self, right_object))

        return YTArray(super(YTArray, self).__add__(right_object))

    def __radd__(self, left_object):
        """ See __add__. """
        if isinstance(left_object, YTArray):
            if not self.units.same_dimensions_as(left_object.units):
                raise UnitOperationError("You cannot add these YTArrays "
                    "because the unit dimensions do not match. `%s + %s` is "
                    "ill-defined." % (left_object.units, self.units))
        else:
            if not self.units.is_dimensionless:
                raise UnitOperationError("You cannot add a pure number to a "
                    "dimensional YTArray. `%s + %s` is ill-defined."
                    % (left_object, self))

        return YTArray(super(YTArray, self).__radd__(left_object))

    def __iadd__(self, other):
        """ See __add__. """
        return np.add(self, other, out=self)

    def __sub__(self, right_object):
        """
        Subtract the object on the right of the `-` from this ytarray. Must
        check for the correct (same dimension) units.

        """
        # Make sure the other object is a YTArray before we use the `units`
        # attribute.
        if isinstance(right_object, YTArray):
            if not self.units.same_dimensions_as(right_object.units):
                raise UnitOperationError("You cannot subtract these YTArrays "
                    "because the unit dimensions do not match. `%s + %s` is "
                    "ill-defined." % (self.units, right_object.units))
        # If the other object is not a YTArray, the only valid case is adding
        # dimensionless things.
        else:
            if not self.units.is_dimensionless:
                raise UnitOperationError("You cannot subtract a pure number "
                    "from a dimensional YTArray. `%s + %s` is ill-defined."
                    % (self, right_object))

        return YTArray(super(YTArray, self).__sub__(right_object))

    def __rsub__(self, left_object):
        """ See __sub__. """
        if isinstance(left_object, YTArray):
            if not self.units.same_dimensions_as(left_object.units):
                raise UnitOperationError("You cannot subtract these YTArrays "
                    "because the unit dimensions do not match. `%s - %s` is "
                    "ill-defined." % (left_object.units, self.units))
        else:
            if not self.units.is_dimensionless:
                raise UnitOperationError("You cannot subtract a dimensional "
                    "from a pure number. `%s - %s` is ill-defined."
                    % (left_object, self))

        return YTArray(super(YTArray, self).__rsub__(left_object))

    def __isub__(self, other):
        """ See __sub__. """
        return np.subtract(self, other, out=self)

    def __neg__(self):
        """ Negate the data. """
        return YTArray(super(YTArray, self).__neg__())

    def __mul__(self, right_object):
        """
        Multiply this YTArray by the object on the right of the `*` operator.
        The unit objects handle being multiplied.

        """
        return YTArray(super(YTArray, self).__mul__(right_object))

    def __rmul__(self, left_object):
        """ See __mul__. """
        return YTArray(super(YTArray, self).__rmul__(left_object))

    def __imul__(self, other):
        """ See __mul__. """
        return np.multiply(self, other, out=self)

    def __div__(self, right_object):
        """
        Divide this YTArray by the object on the right of the `/` operator.

        """
        return YTArray(super(YTArray, self).__div__(right_object))

    def __rdiv__(self, left_object):
        """ See __div__. """
        return YTArray(super(YTArray, self).__rdiv__(left_object))

    def __idiv__(self, other):
        """ See __div__. """
        return np.divide(self, other, out=self)

    def __pow__(self, power):
        """
        Raise this YTArray to some power.

        Parameters
        ----------
        power : float or dimensionless YTArray.
            The pow value.

        """
        if isinstance(power, YTArray):
            if not power.units.is_dimensionless:
                raise UnitOperationError("The power argument must be "
                    "dimensionless. (%s)**(%s) is ill-defined." % (self, power))

        return YTArray(super(YTArray, self).__pow__(power))

    def __abs__(self):
        """ Return a YTArray with the abs of the data. """
        return YTArray(super(YTArray, self).__abs__())

    def sqrt(self):
        """
        Return sqrt of this YTArray. We take the sqrt for the array and use
        take the 1/2 power of the units.

        """
        return YTArray(super(YTArray, self).sqrt(),
                       input_units=self.units**0.5)

    def exp(self):
        """
        Return exp of this YTArray. Checks that the YTArray is dimensionless.

        """
        if not self.units.is_dimensionless:
            raise UnitOperationError("The argument of an exponential must be "
                "dimensionless. exp(%s) is ill-defined." % self)

        return YTArray(super(YTArray, self).exp())

    #
    # Start comparison operators.
    #

    # @todo: outsource to a single method with an op argument.
    def __lt__(self, other):
        """ Test if this is less than the object on the right. """
        # Check that other is a YTArray.
        if isinstance(other, YTArray):
            if not self.units.same_dimensions_as(other.units):
                raise UnitOperationError("The less than operator for "
                    "YTArrays with units %s and %s is not well defined."
                    % (self.units, other.units))

            return super(YTArray, self).__lt__(other.in_units(self.units))

        return super(YTArray, self).__lt__(other)

    def __le__(self, other):
        """ Test if this is less than or equal to the object on the right. """
        # Check that other is a YTArray.
        if isinstance(other, YTArray):
            if not self.units.same_dimensions_as(other.units):
                raise UnitOperationError("The less than or equal operator for "
                    "YTArrays with units %s and %s is not well defined."
                    % (self.units, other.units))

            return super(YTArray, self).__le__(right_object.in_units(self.units))

        return super(YTArray, self).__le__(right_object)

    def __eq__(self, other):
        """ Test if this is equal to the object on the right. """
        # Check that other is a YTArray.
        if isinstance(other, YTArray):
            if not self.units.same_dimensions_as(other.units):
                raise UnitOperationError("The equal operator for "
                    "YTArrays with units %s and %s is not well defined."
                    % (self.units, other.units))

            return super(YTArray, self).__eq__(other.in_units(self.units))

        return super(YTArray, self).__eq__(other)

    def __ne__(self, other):
        """ Test if this is not equal to the object on the right. """
        # Check that the other is a YTArray.
        if isinstance(other, YTArray):
            if not self.units.same_dimensions_as(other.units):
                raise Exception("The not equal operator for "
                    "YTArrays with units %s and %s is not well defined."
                    % (self.units, other.units))

            return super(YTArray, self).__ne__(other.in_units(self.units))

        return super(YTArray, self).__ne__(other)

    def __ge__(self, other):
        """ Test if this is greater than or equal to other. """
        # Check that the other is a YTArray.
        if isinstance(right_object, YTArray):
            if not self.units.same_dimensions_as(other.units):
                raise Exception("The greater than or equal operator for "
                    "YTArrays with units %s and %s is not well defined."
                    % (self.units, other.units))

            return super(YTArray, self).__ge__(other.in_units(self.units))

        return super(YTArray, self).__ge__(other)

    def __gt__(self, other):
        """ Test if this is greater than the object on the right. """
        # Check that the other is a YTArray.
        if isinstance(other, YTArray):
            if not self.units.same_dimensions_as(other.units):
                raise Exception("The greater than operator for "
                    "YTArrays with units %s and %s is not well defined."
                    % (self.units, other.units))

            return super(YTArray, self).__gt__(other.in_units(self.units))

        return super(YTArray, self).__gt__(other)

    #
    # End comparison operators
    #

    def __array_wrap__(self, out_arr, context=None):
        if context is None:
            pass
        elif len(context[1]) == 1:
            # unary operators
            out_arr.units = self._ufunc_registry[context[0]](context[1][0].units)
        elif len(context[1]) == 2:
            # binary operators
            try:
                unit1 = context[1][0].units
            except AttributeError:
                unit1 = Unit()
            try:
                unit2 = context[1][1].units
            except AttributeError:
                if context[0] is power:
                    unit2 = context[1][1]
                else:
                    unit2 = Unit()
            out_arr.units = self._ufunc_registry[context[0]](unit1, unit2)
        else:
            raise RuntimeError("Only unary and binary operators are allowed.")
        if out_arr.size == 1 and out_arr.size != self.size:
            return out_arr[0]
        return super(YTArray, self).__array_wrap__(out_arr, context)
