"""
IndexArray class.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------


import numpy as np

from yt.units.yt_array import \
    binary_operators, \
    commutative_operators, \
    ensure_tuple, \
    iterable, \
    NULL_UNIT, \
    same_unit_operators, \
    UFUNC_REGISTRY, \
    unary_operators, \
    YTArray, \
    YTUfuncUnitError, \
    YTUnitOperationError, \
    YTQuantity
from yt.units.unit_object import \
    Unit
from yt.utilities.exceptions import \
    YTImmutableUnitsError

ELLIPSIS_TYPE = type(Ellipsis)

class IndexArray(YTArray):
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

        if not iterable(input_units):
            input_units = (input_units, )*ndim

        if ndim != len(input_units):
            raise RuntimeError

        obj.units = tuple(Unit(iu, registry=registry) if not isinstance(iu, Unit)
                          else iu for iu in input_units)

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
        self.units = getattr(obj, 'units', (NULL_UNIT, )*ndim)

    def __getitem__(self, item):
        ret = super(IndexArray, self).__getitem__(item)
        if iterable(item):
            if isinstance(item[0], (slice, ELLIPSIS_TYPE)):
                ret = YTArray(ret.view(np.ndarray), self.units[item[1]])
            elif isinstance(item[1], (slice, ELLIPSIS_TYPE)):
                # ret maintains units of original array
                pass
            else:
                ret = YTQuantity(ret, self.units[item[1]])
        else:
            # If we are just getting *one* item back
            if ret.size == 1:
                item_ind = item % self.shape[-1]
                ret = YTQuantity(ret, self.units[item_ind])
            else:
                # this case selects a single row, so we maintain the units
                # of the original array
                pass
        return ret

    def __array_wrap__(self, out_arr, context=None):
        # note that we explicitly call YTArray's superclass, not IndexArray
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
        IndexArray

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

        new_array = IndexArray(new_values, new_units)

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
