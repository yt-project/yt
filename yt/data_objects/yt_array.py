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

from yt.utilities.units import Unit

class UnitOperationError(Exception):
    pass

class YTArray(np.ndarray):
    """

    """
    def __new__(cls, input_array, input_units=None):
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
                raise UnitOperationError("You cannot add these quantities because their dimensions do not match. `%s + %s` is ill-defined" % (self.units, right_object.units))
        else:
            # case of dimensionless self + float
            # the only way this works is with a float so...
            if not self.units.is_dimensionless:
                raise Exception("You cannot add a pure number to a dimensional quantity. `%s + %s` is ill-defined." % (self, right_object))

        # `get_data_in` will not apply the conversion if the units are the same
        return YTArray(super(YTArray, self).__add__(right_object), input_units = right_object.units)

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
        return YTArray(super(YTArray, self).__radd__(left_object), input_units = right_object.units)

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
        return YTArray(super(YTArray, self).__sub__(right_object), input_units = right_object.units)


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
        return YTArray(super(YTArray, self).__rsub__(left_object), input_units = left_object.units)

    def __neg__(self):
        """ Negate the data. """
        return YTArray(super(YTArray, self).__neg__(self), input_units = self.units)

    def __mul__(self, right_object):
        """
        Multiply this quantity by the object on the right of the `*` operator.
        The unit objects handle being multiplied by each other.

        """
        if isinstance(right_object, YTArray):
            return YTArray(super(YTArray, self).__mul__(right_object),
                           input_units=(self.units * right_object.units))

        # `right_object` is not a Quantity object, so try to use it as
        # dimensionless data.
        return YTArray(super(YTArray, self).__mul__(right_object),
                       input_units=self.units)

    def __rmul__(self, left_object):
        """
        Multiply this quantity by the object on the left of the `*` operator.
        The unit objects handle being multiplied by each other.

        """
        if isinstance(left_object, Quantity):
            return Quantity(left_object.data * self.data,
                            left_object.units * self.units)

        # `left_object` is not a Quantity object, so try to use it as
        # dimensionless data.
        return Quantity(left_object * self.data, self.units)

    def __div__(self, right_object):
        """
        Divide this quantity by the object on the right of the `/` operator. The
        unit objects handle being divided by each other.

        """
        if isinstance(right_object, Quantity):
            return Quantity(self.data / right_object.data,
                            self.units / right_object.units)

        # `right_object` is not a Quantity object, so try to use it as
        # dimensionless data.
        return Quantity(self.data / right_object, self.units)

    def __rdiv__(self, left_object):
        """
        Divide the object on the left of the `/` operator by this quantity. The
        unit objects handle being divided by each other.

        """
        if isinstance(left_object, Quantity):
            return Quantity(left_object.data / self.data,
                            left_object.units / self.units)

        # `left_object` is not a Quantity object, so try to use it as
        # dimensionless data.
        return Quantity(left_object / self.data, self.units**(-1))

    def __pow__(self, power):
        """
        Raise this quantity to some power.

        Parameters
        ----------
        power : float or dimensionless Quantity object
            The pow value.

        """
        if isinstance(power, Quantity):
            if power.units.is_dimensionless:
                return Quantity(self.data**power.data, self.units**power.data)
            else:
                raise Exception("The power argument must be dimensionless. (%s)**(%s) is ill-defined." % (self, power))

        return Quantity(self.data**power, self.units**power)

    def __abs__(self):
        """ Return a Quantity with the abs of the data. """
        return Quantity(abs(self.data), self.units)

    def sqrt(self):
        """
        Return sqrt of this Quantity. This is just a wrapper of Quantity.__pow__
        for numpy.sqrt.

        """
        return self**(1.0/2)

    def exp(self):
        """
        Return exp of this Quantity. Ensures that Quantity is dimensionless,
        like __pow__.

        """
        if not self.units.is_dimensionless:
            raise Exception("The argument of an exponential must be dimensionless. exp(%s) is ill-defined." % self)

        try:
            from numpy import exp
        except ImportError:
            raise Exception("This method requires the numpy package. Please install it before calling exp(Quantity)")

        return exp(self.data)

    ### comparison operators
    # @todo: outsource to a single method with an op argument.
    def __lt__(self, right_object):
        """ Test if this is less than the object on the right. """
        # Check that the other is a Quantity.
        if not isinstance(right_object, Quantity):
            raise Exception("You cannot compare a Quantity to a non-Quantity object. %s < %s is ill-defined." % (self, right_object))
        # Check that the dimensions are the same.
        if not self.units.same_dimensions_as(right_object.units):
            raise Exception("You cannot compare quantities of units %s and %s." % (self.units, right_object.units))

        if self.data < right_object.get_data_in(self.units):
            return True
        return False

    def __le__(self, right_object):
        """ Test if this is less than or equal to the object on the right. """
        # Check that the other is a Quantity.
        if not isinstance(right_object, Quantity):
            raise Exception("You cannot compare a Quantity to a non-Quantity object. %s <= %s is ill-defined." % (self, right_object))
        # Check that the dimensions are the same.
        if not self.units.same_dimensions_as(right_object.units):
            raise Exception("You cannot compare quantities of units %s and %s." % (self.units, right_object.units))

        if self.data <= right_object.get_data_in(self.units):
            return True
        return False

    def __eq__(self, right_object):
        """ Test if this is equal to the object on the right. """
        # Check that the other is a Quantity.
        if not isinstance(right_object, Quantity):
            raise Exception("You cannot compare a Quantity to a non-Quantity object. %s == %s is ill-defined." % (self, right_object))
        # Check that the dimensions are the same.
        if not self.units.same_dimensions_as(right_object.units):
            raise Exception("You cannot compare quantities of units %s and %s." % (self.units, right_object.units))

        if self.data == right_object.get_data_in(self.units):
            return True
        return False

    def __ne__(self, right_object):
        """ Test if this is not equal to the object on the right. """
        # Check that the other is a Quantity.
        if not isinstance(right_object, Quantity):
            raise Exception("You cannot compare a Quantity to a non-Quantity object. %s != %s is ill-defined." % (self, right_object))
        # Check that the dimensions are the same.
        if not self.units.same_dimensions_as(right_object.units):
            raise Exception("You cannot compare quantities of units %s and %s." % (self.units, right_object.units))

        if self.data != right_object.get_data_in(self.units):
            return True
        return False

    def __ge__(self, right_object):
        """
        Test if this is greater than or equal to the object on the right.

        """
        # Check that the other is a Quantity.
        if not isinstance(right_object, Quantity):
            raise Exception("You cannot compare a Quantity to a non-Quantity object. %s >= %s is ill-defined." % (self, right_object))
        # Check that the dimensions are the same.
        if not self.units.same_dimensions_as(right_object.units):
            raise Exception("You cannot compare quantities of units %s and %s." % (self.units, right_object.units))

        if self.data >= right_object.get_data_in(self.units):
            return True
        return False

    def __gt__(self, right_object):
        """ Test if this is greater than the object on the right. """
        # Check that the other is a Quantity.
        if not isinstance(right_object, Quantity):
            raise Exception("You cannot compare a Quantity to a non-Quantity object. %s > %s is ill-defined." % (self, right_object))
        # Check that the dimensions are the same.
        if not self.units.same_dimensions_as(right_object.units):
            raise Exception("You cannot compare quantities of units %s and %s." % (self.units, right_object.units))

        if self.data > right_object.get_data_in(self.units):
            return True
        return False
