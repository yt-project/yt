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


class IndexArray(np.ndarray):
    """
    An ndarray subclass useful for indexing along dimensions that have units
    """

    def __new__(cls, input_array, input_units=None, registry=None, dtype=None):
        if input_array is NotImplemented:
            return input_array
        if isinstance(input_array, YTArray):
            raise NotImplementedError
        elif iterable(input_array):
            pass
        else:
            raise NotImplementedError

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
            if isinstance(item[0], slice):
                ret = YTArray(ret, self.units[item[1]])
            elif isinstance(item[1], slice):
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
        ret = super(IndexArray, self).__array_wrap__(out_arr, context)
        if context is None:
            return ret
        elif context[0] in unary_operators:
            u = context[1][0].units
            units = UFUNC_REGISTRY[context[0]](u)
            if units is None:
                units = (NULL_UNIT, )*len(u)
            ret_class = type(self)
        elif context[0] in binary_operators:
            oper1 = context[1][0]
            oper2 = context[1][1]
            unit1 = getattr(oper1, 'units', None)
            unit2 = getattr(oper2, 'units', None)

            # eventually implement something like get_binary_op_return_class
            ret_class = type(self)

            if unit1 is None:
                unit1 = NULL_UNIT
            if unit2 is None:
                unit2 = NULL_UNIT
            if not iterable(unit1):
                unit1 = (unit1,) * len(unit2)
            if not iterable(unit2):
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

        return ret_class(np.array(out_arr, copy=False), units)

    def __eq__(self, other):
        """ Test if this is equal to the object on the right. """
        if other is None:
            return False
        if not hasattr(other, 'units'):
            if not all([u == NULL_UNIT for u in self.units]):
                return False
        if not iterable(other.units):
            return False
        if not all([u1 == u2 for u1, u2 in zip(self.units, other.units)]):
            return False
        return np.array(self).__eq__(np.array(other))

    def __str__(self):
        return super(IndexArray, self).__str__() + ' ' + str(self.units)

    def __repr__(self):
        return super(IndexArray, self).__repr__() + ' ' + repr(self.units)

    def ndarray_view(self):
        """
        Returns a view into the array, but as an ndarray rather than ytarray.

        Returns
        -------
        View of this array's data.
        """
        return self.view(np.ndarray)

    @property
    def ndview(self):
        """Get a view of the array data."""
        return self.ndarray_view()

    d = ndview
