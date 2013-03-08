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
    negative, absolute, sqrt, square, power, reciprocal

from yt.utilities.units import Unit, dimensionless

class UnitOperationError(Exception):
    pass

def sqrt_unit(unit):
    return unit**0.5

def multiply_units(unit1, unit2):
    return unit1 * unit2

def preserve_unit(unit1, unit2):
    return unit1

def power_unit(unit, power):
    return unit**power

def divide_units(unit1, unit2):
    return unit1/unit2

class YTArray(np.ndarray):
    """

    """
    _ufunc_registry = {sqrt: sqrt_unit, multiply: multiply_units,
                       add: preserve_unit, subtract: preserve_unit,
                       power: power_unit, divide: divide_units}

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
        if not isinstance(units, Unit):
            units = Unit(units)

        if not self.units.same_dimensions_as(units):
            raise Exception("Cannot convert to units with different dimensionality. Current unit is %s, argument is %s" % (self.units.dimensions, units))

        return units

    def convert_to(self, units):
        """
        Convert the data and units to given unit. This overwrites the ``data``
        and ``units`` attributes, making no copies, and returns None.

        Parameters
        ----------
        units : Unit object or string
            The units you want the data in.

        """
        new_units = self._unit_repr_check_same(units)
        conversion_factor = get_conversion_factor(self.units, new_units)
        self.data *= conversion_factor
        self.units = new_units

        return self

    def convert_to_cgs(self):
        """
        Convert the data and units to the equivalent cgs units. This overwrites
        the ``data`` and ``units`` attributes, making no copies, and returns
        None.

        """
        return self.convert_to(self.units.get_cgs_equivalent())

    def get_in(self, units):
        """
        Creates a new Quantity with the data in the supplied units, and returns
        it. Does not modify this object.

        Parameters
        ----------
        units : Unit object or string
            The units you want to get a new quantity in.

        Returns
        -------
        Quantity object with converted data and supplied units.

        """
        new_units = self._unit_repr_check_same(units)
        conversion_factor = get_conversion_factor(self.units, new_units)

        return Quantity(self.data * conversion_factor, new_units)

    def get_in_cgs(self):
        """
        Creates a new Quantity with the data in the equivalent cgs units, and
        returns it. Does not modify this object.

        Returns
        -------
        Quantity object with data converted to cgs and cgs units.

        """
        return self.get_in(self.units.get_cgs_equivalent())

    def get_data_in(self, units):
        """
        Returns the data, converted to the supplied units.

        Parameters
        ----------
        units : Unit object or string
            The units you want the data in.

        Returns
        -------
        ``data`` attribute, multiplied by the conversion factor to the supplied
        units.

        """
        new_units = self._unit_repr_check_same(units)

        # don't operate on data if given the same units
        if self.units == new_units:
            return self.data

        conversion_factor = get_conversion_factor(self.units, new_units)

        return self.data * conversion_factor

    def get_data_in_cgs(self):
        """
        Returns the data, multiplied by the conversion factor to cgs.

        """
        return self.get_data_in(self.units.get_cgs_equivalent())

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
        if isinstance(right_object, YTArray):
            # make sure it's a quantity before we check units attribute
            if not self.units.same_dimensions_as(right_object.units):
                raise UnitOperationError("You cannot add these quantities because "
                                         "their dimensions do not match. "
                                         "`%s + %s` is ill-defined" % (self.units, right_object.units))
        else:
            # case of dimensionless self + float
            # the only way this works is with a float so...
            if not self.units.is_dimensionless:
                raise Exception("You cannot add a pure number to a "
                                "dimensional quantity. `%s + %s` is ill-defined." % (self, right_object))

        # `get_data_in` will not apply the conversion if the units are the same
        return YTArray(super(YTArray, self).__add__(right_object))

    def __radd__(self, left_object):
        """
        Add this ytarray to the object on the left of the `+` operator. Must
        check for the correct (same dimension) units.

        """
        if isinstance(left_object, YTArray):
            # make sure it's a quantity before we check units attribute
            if not self.units.same_dimensions_as(left_object.units):
                raise Exception("You cannot add these quantities because their dimensions do not match. `%s + %s` is ill-defined" % (left_object.units, self.units))
        else:  # the only way this works is with a float so...
            # case of dimensionless float + self
            if not self.units.is_dimensionless:
                raise Exception("You cannot add a pure number to a dimensional quantity. `%s + %s` is ill-defined." % (left_object, self))

        # `get_data_in` will not apply the conversion if the units are the same
        return YTArray(super(YTArray, self).__radd__(left_object))

    def __iadd__(self, other):
        self = self + other
        return self

    def __sub__(self, right_object):
        """
        Subtract the object on the right of the `-` from this quantity. Must
        check for the correct (same dimension) units. If the quantities have
        different units, we always use the units on the left.

        """
        if isinstance(right_object, YTArray):  # make sure it's a quantity before we check units attribute
            if not self.units.same_dimensions_as(right_object.units):
                raise Exception("You cannot add these quantities because their dimensions do not match. `%s - %s` is ill-defined" % (self.units, right_object.units))
        else:
            # case of dimensionless self + float
            # the only way this works is with a float so...
            if not self.units.is_dimensionless:
                raise Exception("You cannot add a pure number to a dimensional quantity. `%s - %s` is ill-defined." % (self, right_object))

        # `get_data_in` will not apply the conversion if the units are the same
        return YTArray(super(YTArray, self).__sub__(right_object))
        
    def __rsub__(self, left_object):
        """
        Subtract this quantity from the object on the left of the `-` operator.
        Must check for the correct (same dimension) units. If the quantities
        have different units, we always use the units on the left.

        """
        if isinstance(left_object, Quantity):  # make sure it's a quantity before we check units attribute
            if not self.units.same_dimensions_as(left_object.units):
                raise Exception("You cannot add these quantities because their dimensions do not match. `%s - %s` is ill-defined" % (left_object.units, self.units))
        else:
            # case of dimensionless self + float
            # the only way this works is with a float so...
            if not self.units.is_dimensionless:
                raise Exception("You cannot add a pure number to a dimensional quantity. `%s - %s` is ill-defined." % (left_object, self))

        # `get_data_in` will not apply the conversion if the units are the same
        return YTArray(super(YTArray, self).__rsub__(left_object))

    def __isub__(self, other):
        self = self - other
        return self

    def __neg__(self):
        """ Negate the data. """
        return YTArray(super(YTArray, self).__neg__(self))

    def __mul__(self, right_object):
        """
        Multiply this quantity by the object on the right of the `*` operator.
        The unit objects handle being multiplied by each other.

        """
        if isinstance(right_object, YTArray):
            return YTArray(super(YTArray, self).__mul__(right_object))

        # `right_object` is not a Quantity object, so try to use it as
        # dimensionless data.
        return YTArray(super(YTArray, self).__mul__(right_object))

    def __rmul__(self, left_object):
        """
        Multiply this quantity by the object on the left of the `*` operator.
        The unit objects handle being multiplied by each other.

        """
        if isinstance(left_object, YTArray):
            return YTArray(super(YTArray, self).__rmul__(left_object))

        # `left_object` is not a Quantity object, so try to use it as
        # dimensionless data.
        return YTArray(super(YTArray, self).__rmul__(left_object))

    def __imul__(self, other):
        self = self * other
        return self

    def __div__(self, right_object):
        """
        Divide this quantity by the object on the right of the `/` operator. The
        unit objects handle being divided by each other.

        """
        if isinstance(right_object, YTArray):
            return YTArray(super(YTArray, self).__div__(right_object))

        # `right_object` is not a  object, so try to use it as
        # dimensionless data.
        return YTArray(super(YTArray, self).__div__(right_object))

    def __rdiv__(self, left_object):
        """
        Divide the object on the left of the `/` operator by this quantity. The
        unit objects handle being divided by each other.

        """
        if isinstance(left_object, YTArray):
            return YTArray(super(YTArray, self).__rdiv__(left_object))

        # `left_object` is not a Quantity object, so try to use it as
        # dimensionless data.
        return YTArray(super(YTArray, self).__rdiv__(left_object))

    def __idiv__(self, other):
        self = self / other
        return self

    def __pow__(self, power):
        """
        Raise this quantity to some power.

        Parameters
        ----------
        power : float or dimensionless Quantity object
            The pow value.

        """
        if isinstance(power, YTArray):
            if not power.units.is_dimensionless:
                raise Exception("The power argument must be dimensionless. (%s)**(%s) is ill-defined." % (self, power))

        return YTArray(super(YTArray, self).__pow__(power))

    def __abs__(self):
        """ Return a Quantity with the abs of the data. """
        return YTArray(super(YTArray, self).__abs__())
                 
    def sqrt(self):
        """
        Return sqrt of this Quantity. This is just a wrapper of Quantity.__pow__
        for numpy.sqrt.

        """
        return YTArray(super(YTArray, self).sqrt(), input_units = self.units**0.5)

    def exp(self):
        """
        Return exp of this Quantity. Ensures that Quantity is dimensionless,
        like __pow__.

        """
        if not self.units.is_dimensionless:
            raise Exception("The argument of an exponential must be dimensionless. exp(%s) is ill-defined." % self)

        return exp(self.data)

    ### comparison operators
    # @todo: outsource to a single method with an op argument.
    def __lt__(self, right_object):
        """ Test if this is less than the object on the right. """
        # Check that the other is a YTArray.
        if isinstance(right_object, YTArray):
            if not self.units.same_dimensions_as(right_object.units):
                raise Exception("The less than operator for quantities "
                                "with units %s and %s is not well defined." 
                                % (self.units, right_object.units))

            return super(YTArray, self).__lt__(right_object.get_data_in(self.units))
        return super(YTArray, self).__lt__(right_object)

    def __le__(self, right_object):
        """ Test if this is less than or equal to the object on the right. """
        # Check that the other is a YTArray.
        if isinstance(right_object, Quantity):
            if not self.units.same_dimensions_as(right_object.units):
                raise Exception("The less than or equal operator for quantities "
                                "with units %s and %s is not well defined." 
                                % (self.units, right_object.units))

            return super(YTArray, self).__le__(right_object.get_data_in(self.units))
        return super(YTArray, self).__lt__(right_object)

    def __eq__(self, right_object):
        """ Test if this is equal to the object on the right. """
        # Check that the other is a YTArray.
        if isinstance(right_object, YTArray):
            if not self.units.same_dimensions_as(right_object.units):
                raise Exception("The equality operator for quantities with units" 
                                "%s and %s is not well defined." 
                                % (self.units, right_object.units))

            return super(YTArray, self).__eq__(right_object.get_data_in(self.units))
        return super(YTArray, self).__eq__(right_object)

    def __ne__(self, right_object):
        """ Test if this is not equal to the object on the right. """
        # Check that the other is a YTArray.
        if isinstance(right_object, YTArray):
            if not self.units.same_dimensions_as(right_object.units):
                raise Exception("The not equal operator for quantities "
                                "with units %s and %s is not well defined." 
                                % (self.units, right_object.units))
            
            return super(YTArray, self).__ne__(right_object.get_data_in(self.units))
        return super(YTArray, self).__ne__(right_object)

    def __ge__(self, right_object):
        """
        Test if this is greater than or equal to the object on the right.

        """
        # Check that the other is a YTArray.
        if isinstance(right_object, YTArray):
            if not self.units.same_dimensions_as(right_object.units):
                raise Exception("The greater than or equal operator for quantities "
                                "with units %s and %s is not well defined." 
                                % (self.units, right_object.units))

            return super(YTArray, self).__ge__(right_object.get_data_in(self.units))
        return super(YTArray, self).__ge__(right_objects)

    def __gt__(self, right_object):
        """ Test if this is greater than the object on the right. """
        # Check that the other is a YTArray.
        if isinstance(right_object, YTArray):
            if not self.units.same_dimensions_as(right_object.units):
                raise Exception("The greater than operator for quantities "
                                "with units %s and %s is not well defined." 
                                % (self.units, right_object.units))

            return super(YTArray, self).__gt__(right_object.get_data_in(self.units))
        return super(YTArray, self).__gt__(right_object)

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
        return out_arr
