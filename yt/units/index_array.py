"""
YTIndexArray class.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------


import numpy as np
import operator

from numbers import Number as numeric_type
from six import string_types
from six.moves import reduce

from yt.units.yt_array import \
    binary_operators, \
    commutative_operators, \
    ensure_tuple, \
    iterable, \
    NULL_UNIT, \
    return_arr, \
    same_unit_operators, \
    UFUNC_REGISTRY, \
    unary_operators, \
    YTArray, \
    YTUfuncUnitError, \
    YTUnitOperationError, \
    YTQuantity
from yt.units.unit_object import \
    Unit, UnitTuple
from yt.utilities.exceptions import \
    YTImmutableUnitsError

ELLIPSIS_TYPE = type(Ellipsis)

class YTIndexArray(YTArray):
    """
    An ndarray subclass useful for indexing along dimensions that have units
    """

    __array_priority__ = 3.0

    def __new__(cls, input_array, input_units=None, registry=None,
                dtype=np.float64):
        if input_array is NotImplemented:
            return input_array
        if input_units is None:
            input_units = getattr(input_array, 'units', '')

        # Create an instance of our array subclass
        obj = np.asarray(input_array, dtype=dtype).view(cls)

        if len(obj.shape) == 1:
            ndim = obj.shape[0]
        elif len(obj.shape) == 2:
            ndim = obj.shape[-1]
        else:
            raise NotImplementedError

        if not iterable(input_units) or isinstance(input_units, string_types):
            input_units = (input_units, )*ndim

        if ndim != len(input_units):
            raise RuntimeError(
                'Cannot create index array with shape %s and units %s'
                % (obj.shape, input_units)
            )

        obj.units = UnitTuple(
            (Unit(iu, registry=registry) if not isinstance(iu, Unit)
             else iu for iu in input_units))

        return obj

    def __array_finalize__(self, obj):
        if obj is None and hasattr(self, 'unit'):
            return
        if len(obj.shape) == 1:
            ndim = obj.shape[0]
        elif len(obj.shape) == 2:
            ndim = obj.shape[-1]
        else:
            raise NotImplementedError
        units = getattr(obj, 'units', (NULL_UNIT, )*ndim)
        self.units = units

    def __getslice__(self, i, j):
        # This is only called under python2. It is needed because py2 extension
        # types (like ndarray) still implicitly call __getslice__
        return self[slice(i, j)]

    def __getitem__(self, item):
        ret = self.view(np.ndarray)[item]
        item_type = type(item)
        if ret.size == 1:
            if item_type is slice:
                item_ind = item.start
                item_ind = item_ind % self.shape[-1]
            elif item_type is tuple:
                item_ind = item[1]
            else:
                item_ind = item % self.shape[-1]
            ret_class = YTQuantity
            ret_units = self.units[item_ind]
        else:
            if item_type is ELLIPSIS_TYPE:
                ret_class = YTIndexArray
                ret_units = self.units
            elif item_type is tuple:
                if isinstance(item[0], numeric_type):
                    ret_class = YTIndexArray
                else:
                    if len(ret.shape) == 1 or ret.shape[1] == 1:
                        ret_class = YTArray
                    else:
                        ret_class = YTIndexArray
                ret_units = self.units[item[1]]
            elif item_type is list or issubclass(item_type, np.ndarray):
                ret_class = YTIndexArray
                if len(self.shape) == 2:
                    ret_units = self.units
                else:
                    ret_units = UnitTuple([self.units[i] for i in item])
            else:
                ret_class = YTIndexArray
                if len(self.shape) == 1 and item_type is slice:
                    ret_units = self.units[item]
                else:
                    ret_units = self.units
        return ret_class(ret, input_units=ret_units)

    def __array_wrap__(self, out_arr, context=None):
        # note that we explicitly call YTArray's superclass, not YTIndexArray
        ret = super(YTArray, self).__array_wrap__(out_arr, context)
        if context is None:
            return ret
        elif context[0] in unary_operators:
            u = context[1][0].units
            units = UFUNC_REGISTRY[context[0]](u)
            ret_class = type(self)
        elif context[0] in binary_operators:
            oper1 = context[1][0]
            oper2 = context[1][1]
            unit1 = getattr(oper1, 'units', None)
            unit2 = getattr(oper2, 'units', None)

            # eventually implement something like get_binary_op_return_class
            ret_class = type(self)

            if unit1 is None:
                unit1 = Unit(registry=getattr(unit2[0], 'registry', None))
            elif unit2 is None:
                unit2 = Unit(registry=getattr(unit1[0], 'registry', None))

            if not iterable(unit1):
                unit1 = (unit1,) * len(unit2)
            elif not iterable(unit2):
                unit2 = (unit2,) * len(unit1)

            unit_operator = UFUNC_REGISTRY[context[0]]
            if unit_operator in same_unit_operators:
                if unit1 != unit2:
                    if not unit1.same_dimensions_as(unit2):
                        raise YTUnitOperationError(context[0], unit1, unit2)
                    else:
                        raise YTUfuncUnitError(context[0], unit1, unit2)
            units = tuple(unit_operator(u1, u2) for u1, u2 in zip(unit1, unit2))
            if all([u is None for u in units]):
                units = None
            if unit_operator in commutative_operators:
                for unit in units:
                    if unit.is_dimensionless and unit.base_value != 1.0:
                        # fix this later
                        raise NotImplementedError
        # Deal with units for functions like isnan, which should always return
        # ndarray instances
        if units is None:
            out_arr = np.array(out_arr, copy=False)
            return out_arr
        return ret_class(np.array(out_arr, copy=False), units)

    #
    # Begin unit conversion methods
    #

    def convert_to_units(self, units):
        """
        This function raises a RuntimeError, use in_units instead.

        """
        raise YTImmutableUnitsError('in_units')

    def convert_to_base(self):
        """
        This function raises a RuntimeError, use in_base instead.

        """
        raise YTImmutableUnitsError('in_base')

    def convert_to_cgs(self):
        """
        This function raises a RuntimeError, use in_cgs instead.

        """
        raise YTImmutableUnitsError('in_cgs')

    def convert_to_mks(self):
        """
        This function raises a RuntimeError, use in_mks instead.

        """
        YTImmutableUnitsError('in_mks')

    def in_units(self, units):
        """
        Retuns a copy of this array with the data in the supplied units.

        Parameters
        ----------
        units : Unit object or string
            The units you want to get a new quantity in.

        Returns
        -------
        YTIndexArray

        """
        units = ensure_tuple(units)
        if len(units) == 1:
            units = units * len(self.units)

        new_units = []
        new_values = self.v

        for i, (sunit, unit) in enumerate(zip(self.units, units)):
            nu = sunit._unit_repr_check_same(unit)
            new_units.append(nu)
            conversion_factor, offset = sunit.get_conversion_factor(nu)
            new_values[..., i] *= conversion_factor

        new_array = YTIndexArray(new_values, new_units)

        if offset:
            # punt on this for now
            raise NotImplementedError

        return new_array

    def in_base(self):
        """
        Creates a copy of this array with the data in the equivalent base units,
        and returns it.

        Returns
        -------
        Quantity object with data converted to cgs units.

        """
        return self.in_units(
            tuple([u.get_base_equivalent() for u in self.units])
        )

    def in_cgs(self):
        """
        Creates a copy of this array with the data in the equivalent cgs units,
        and returns it.

        Returns
        -------
        Quantity object with data converted to cgs units.

        """
        return self.in_units(
            tuple([u.get_cgs_equivalent() for u in self.units])
        )

    def in_mks(self):
        """
        Creates a copy of this array with the data in the equivalent mks units,
        and returns it.

        Returns
        -------
        Quantity object with data converted to mks units.

        """
        return self.in_units(
            tuple([u.get_mks_equivalent() for u in self.units])
        )

    #
    # End unit conversion methods
    #

    #
    # Begin helper methods
    #

    @classmethod
    def from_astropy(cls, *args, **kwargs):
        raise NotImplementedError(
            "YTIndexArray does not support conversions to and from astropy"
        )

    def to_astropy(self, *args, **kwargs):
        raise NotImplementedError(
            "YTIndexArray does not support conversions to and from astropy"
        )

    @classmethod
    def from_pint(cls, *args, **kwargs):
        raise NotImplementedError(
            "YTIndexArray does not support conversions to and from pint"
        )

    def to_pint(self, *args, **kwargs):
        raise NotImplementedError(
            "YTIndexArray does not support conversions to and from pint"
        )

    def write_hdf5(self, *args, **kwargs):
        raise NotImplementedError(
            "YTIndexArray does not yet support native HDF5 I/O"
        )

    @classmethod
    def from_hdf5(cls, *args, **kwargs):
        raise NotImplementedError(
            "YTIndexArray does not yet support native HDF5 I/O"
        )

    @property
    def unit_quantity(self):
        if self.units.is_homogenous:
            return YTQuantity(1.0, self.units[0])
        msg = "Cannot convert IndexArray with units '%s' to a quantity"
        raise RuntimeError(msg % (self.units, ))

    uq = unit_quantity

    @return_arr
    def prod(self, axis=None, dtype=None, out=None):
        if len(self.shape) == 2:
            nrows = self.shape[0]
        else:
            nrows = 1
        if axis is None:
            units = reduce(operator.mul, self.units)**nrows
        elif axis == 0:
            units = self.units**nrows
        elif axis == 1:
            units = reduce(operator.mul, self.units)
        return np.asarray(self).prod(axis=axis, dtype=dtype, out=out), units

    def _get_reduction_units_from_axis(self, axis, method):
        if axis in (None, 1):
            if self.units.is_homogenous:
                units = self.units[0]
            else:
                raise RuntimeError(
                    'Cannot use the %s method with "axis=%s" for arrays '
                    'with units %s' % (method, axis, self.units, )
                )
        elif axis == 0:
            units = self.units
        return units

    @return_arr
    def sum(self, axis=None, dtype=None, out=None):
        units = self._get_reduction_units_from_axis(axis, 'sum')
        return np.asarray(self).sum(axis=axis, dtype=dtype, out=out), units

    @return_arr
    def min(self, axis=None, out=None, keepdims=False):
        units = self._get_reduction_units_from_axis(axis, 'min')
        return np.asarray(self).min(axis=axis, out=out, keepdims=keepdims), units

    def argmin(self, axis=None, out=None):
        self._get_reduction_units_from_axis(axis, 'argmin')
        return np.asarray(self).argmin(axis=axis, out=out)

    @return_arr
    def max(self, axis=None, out=None, keepdims=False):
        units = self._get_reduction_units_from_axis(axis, 'max')
        return np.asarray(self).max(axis=axis, out=out, keepdims=keepdims), units

    def argmax(self, axis=None, out=None):
        self._get_reduction_units_from_axis(axis, 'argmax')
        return np.asarray(self).argmax(axis=axis, out=out)

    @return_arr
    def mean(self, axis=None, dtype=None, out=None, keepdims=False):
        units = self._get_reduction_units_from_axis(axis, 'max')
        return np.asarray(self).mean(
            axis=axis, dtype=dtype, out=out, keepdims=keepdims), units

    @return_arr
    def std(self, axis=None, dtype=None, out=None, ddof=0, keepdims=False):
        units = self._get_reduction_units_from_axis(axis, 'std')
        return np.asarray(self).std(
            axis=axis, dtype=dtype, out=out, ddof=ddof, keepdims=keepdims), units
