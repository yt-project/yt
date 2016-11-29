"""
YTArray class.



"""
from __future__ import print_function
#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import copy
import numpy as np

from functools import wraps
from numpy import \
    add, subtract, multiply, divide, logaddexp, logaddexp2, true_divide, \
    floor_divide, negative, power, remainder, mod, absolute, rint, \
    sign, conj, exp, exp2, log, log2, log10, expm1, log1p, sqrt, square, \
    reciprocal, ones_like, sin, cos, tan, arcsin, arccos, arctan, arctan2, \
    hypot, sinh, cosh, tanh, arcsinh, arccosh, arctanh, deg2rad, rad2deg, \
    bitwise_and, bitwise_or, bitwise_xor, invert, left_shift, right_shift, \
    greater, greater_equal, less, less_equal, not_equal, equal, logical_and, \
    logical_or, logical_xor, logical_not, maximum, minimum, fmax, fmin, \
    isreal, iscomplex, isfinite, isinf, isnan, signbit, copysign, nextafter, \
    modf, ldexp, frexp, fmod, floor, ceil, trunc, fabs, spacing

from yt.units.unit_object import Unit, UnitParseError
from yt.units.unit_registry import UnitRegistry
from yt.units.dimensions import \
    angle, \
    current_mks, \
    dimensionless, \
    em_dimensions
from yt.utilities.exceptions import \
    YTUnitOperationError, YTUnitConversionError, \
    YTUfuncUnitError, YTIterableUnitCoercionError, \
    YTInvalidUnitEquivalence, YTEquivalentDimsError
from yt.utilities.lru_cache import lru_cache
from numbers import Number as numeric_type
from yt.utilities.on_demand_imports import _astropy
from sympy import Rational
from yt.units.unit_lookup_table import \
    default_unit_symbol_lut
from yt.units.equivalencies import equivalence_registry
from yt.utilities.logger import ytLogger as mylog
from .pint_conversions import convert_pint_units

NULL_UNIT = Unit()

# redefine this here to avoid a circular import from yt.funcs
def iterable(obj):
    try: len(obj)
    except: return False
    return True

def return_arr(func):
    @wraps(func)
    def wrapped(*args, **kwargs):
        ret, units = func(*args, **kwargs)
        if ret.shape == ():
            return YTQuantity(ret, units)
        else:
            # This could be a subclass, so don't call YTArray directly.
            return type(args[0])(ret, units)
    return wrapped

@lru_cache(maxsize=128, typed=False)
def sqrt_unit(unit):
    return unit**0.5

@lru_cache(maxsize=128, typed=False)
def multiply_units(unit1, unit2):
    return unit1 * unit2

def preserve_units(unit1, unit2):
    return unit1

@lru_cache(maxsize=128, typed=False)
def power_unit(unit, power):
    return unit**power

@lru_cache(maxsize=128, typed=False)
def square_unit(unit):
    return unit*unit

@lru_cache(maxsize=128, typed=False)
def divide_units(unit1, unit2):
    return unit1/unit2

@lru_cache(maxsize=128, typed=False)
def reciprocal_unit(unit):
    return unit**-1

def passthrough_unit(unit):
    return unit

def return_without_unit(unit):
    return None

def arctan2_unit(unit1, unit2):
    return NULL_UNIT

def comparison_unit(unit1, unit2):
    return None

def invert_units(unit):
    raise TypeError(
        "Bit-twiddling operators are not defined for YTArray instances")

def bitop_units(unit1, unit2):
    raise TypeError(
        "Bit-twiddling operators are not defined for YTArray instances")

def coerce_iterable_units(input_object):
    if isinstance(input_object, np.ndarray):
        return input_object
    if iterable(input_object):
        if any([isinstance(o, YTArray) for o in input_object]):
            ff = getattr(input_object[0], 'units', NULL_UNIT, )
            if any([ff != getattr(_, 'units', NULL_UNIT) for _ in input_object]):
                raise YTIterableUnitCoercionError(input_object)
            # This will create a copy of the data in the iterable.
            return YTArray(input_object)
        return input_object
    else:
        return input_object

def sanitize_units_mul(this_object, other_object):
    inp = coerce_iterable_units(this_object)
    ret = coerce_iterable_units(other_object)
    # If the other object is a YTArray and has the same dimensions as the object
    # under consideration, convert so we don't mix units with the same
    # dimensions.
    if isinstance(ret, YTArray):
        if inp.units.same_dimensions_as(ret.units):
            ret.in_units(inp.units)
    return ret

def sanitize_units_add(this_object, other_object, op_string):
    inp = coerce_iterable_units(this_object)
    ret = coerce_iterable_units(other_object)
    # Make sure the other object is a YTArray before we use the `units`
    # attribute.
    if isinstance(ret, YTArray):
        if not inp.units.same_dimensions_as(ret.units):
            # handle special case of adding or subtracting with zero or
            # array filled with zero
            if not np.any(other_object):
                return ret.view(np.ndarray)
            elif not np.any(this_object):
                return ret
            raise YTUnitOperationError(op_string, inp.units, ret.units)
        ret = ret.in_units(inp.units)
    else:
        # If the other object is not a YTArray, then one of the arrays must be
        # dimensionless or filled with zeros
        if not inp.units.is_dimensionless and np.any(ret):
            raise YTUnitOperationError(op_string, inp.units, dimensionless)
    return ret

def validate_comparison_units(this, other, op_string):
    # Check that other is a YTArray.
    if hasattr(other, 'units'):
        if this.units.expr is other.units.expr:
            if this.units.base_value == other.units.base_value:
                return other
        if not this.units.same_dimensions_as(other.units):
            raise YTUnitOperationError(op_string, this.units, other.units)
        return other.in_units(this.units)

    return other

@lru_cache(maxsize=128, typed=False)
def _unit_repr_check_same(my_units, other_units):
    """
    Takes a Unit object, or string of known unit symbol, and check that it
    is compatible with this quantity. Returns Unit object.

    """
    # let Unit() handle units arg if it's not already a Unit obj.
    if not isinstance(other_units, Unit):
        other_units = Unit(other_units, registry=my_units.registry)

    equiv_dims = em_dimensions.get(my_units.dimensions, None)
    if equiv_dims == other_units.dimensions:
        if current_mks in equiv_dims.free_symbols:
            base = "SI"
        else:
            base = "CGS"
        raise YTEquivalentDimsError(my_units, other_units, base)

    if not my_units.same_dimensions_as(other_units):
        raise YTUnitConversionError(
            my_units, my_units.dimensions, other_units, other_units.dimensions)

    return other_units

unary_operators = (
    negative, absolute, rint, ones_like, sign, conj, exp, exp2, log, log2,
    log10, expm1, log1p, sqrt, square, reciprocal, sin, cos, tan, arcsin,
    arccos, arctan, sinh, cosh, tanh, arcsinh, arccosh, arctanh, deg2rad,
    rad2deg, invert, logical_not, isreal, iscomplex, isfinite, isinf, isnan,
    signbit, floor, ceil, trunc, modf, frexp, fabs, spacing
)

binary_operators = (
    add, subtract, multiply, divide, logaddexp, logaddexp2, true_divide, power,
    remainder, mod, arctan2, hypot, bitwise_and, bitwise_or, bitwise_xor,
    left_shift, right_shift, greater, greater_equal, less, less_equal,
    not_equal, equal, logical_and, logical_or, logical_xor, maximum, minimum,
    fmax, fmin, copysign, nextafter, ldexp, fmod,
)

trigonometric_operators = (
    sin, cos, tan,
)

class YTArray(np.ndarray):
    """
    An ndarray subclass that attaches a symbolic unit object to the array data.

    Parameters
    ----------

    input_array : iterable
        A tuple, list, or array to attach units to
    input_units : String unit specification, unit symbol object, or astropy units
        The units of the array. Powers must be specified using python
        syntax (cm**3, not cm^3).
    registry : A UnitRegistry object
        The registry to create units from. If input_units is already associated
        with a unit registry and this is specified, this will be used instead of
        the registry associated with the unit object.
    dtype : string or NumPy dtype object
        The dtype of the array data. Defaults to the dtype of the input data,
        or, if none is found, uses np.float64
    bypass_validation : boolean
        If True, all input validation is skipped. Using this option may produce
        corrupted, invalid units or array data, but can lead to significant
        speedups in the input validation logic adds significant overhead. If set,
        input_units *must* be a valid unit object. Defaults to False.

    Examples
    --------

    >>> from yt import YTArray
    >>> a = YTArray([1, 2, 3], 'cm')
    >>> b = YTArray([4, 5, 6], 'm')
    >>> a + b
    YTArray([ 401.,  502.,  603.]) cm
    >>> b + a
    YTArray([ 4.01,  5.02,  6.03]) m

    NumPy ufuncs will pass through units where appropriate.

    >>> import numpy as np
    >>> a = YTArray(np.arange(8), 'g/cm**3')
    >>> np.ones_like(a)
    YTArray([1, 1, 1, 1, 1, 1, 1, 1]) g/cm**3

    and strip them when it would be annoying to deal with them.

    >>> np.log10(a)
    array([       -inf,  0.        ,  0.30103   ,  0.47712125,  0.60205999,
            0.69897   ,  0.77815125,  0.84509804])

    YTArray is tightly integrated with yt datasets:

    >>> import yt
    >>> ds = yt.load('IsolatedGalaxy/galaxy0030/galaxy0030')
    >>> a = ds.arr(np.ones(5), 'code_length')
    >>> a.in_cgs()
    YTArray([  3.08600000e+24,   3.08600000e+24,   3.08600000e+24,
             3.08600000e+24,   3.08600000e+24]) cm

    This is equivalent to:

    >>> b = YTArray(np.ones(5), 'code_length', registry=ds.unit_registry)
    >>> np.all(a == b)
    True

    """
    _ufunc_registry = {
        add: preserve_units,
        subtract: preserve_units,
        multiply: multiply_units,
        divide: divide_units,
        logaddexp: return_without_unit,
        logaddexp2: return_without_unit,
        true_divide: divide_units,
        floor_divide: divide_units,
        negative: passthrough_unit,
        power: power_unit,
        remainder: preserve_units,
        mod: preserve_units,
        fmod: preserve_units,
        absolute: passthrough_unit,
        fabs: passthrough_unit,
        rint: return_without_unit,
        sign: return_without_unit,
        conj: passthrough_unit,
        exp: return_without_unit,
        exp2: return_without_unit,
        log: return_without_unit,
        log2: return_without_unit,
        log10: return_without_unit,
        expm1: return_without_unit,
        log1p: return_without_unit,
        sqrt: sqrt_unit,
        square: square_unit,
        reciprocal: reciprocal_unit,
        ones_like: passthrough_unit,
        sin: return_without_unit,
        cos: return_without_unit,
        tan: return_without_unit,
        sinh: return_without_unit,
        cosh: return_without_unit,
        tanh: return_without_unit,
        arcsin: return_without_unit,
        arccos: return_without_unit,
        arctan: return_without_unit,
        arctan2: arctan2_unit,
        arcsinh: return_without_unit,
        arccosh: return_without_unit,
        arctanh: return_without_unit,
        hypot: preserve_units,
        deg2rad: return_without_unit,
        rad2deg: return_without_unit,
        bitwise_and: bitop_units,
        bitwise_or: bitop_units,
        bitwise_xor: bitop_units,
        invert: invert_units,
        left_shift: bitop_units,
        right_shift: bitop_units,
        greater: comparison_unit,
        greater_equal: comparison_unit,
        less: comparison_unit,
        less_equal: comparison_unit,
        not_equal: comparison_unit,
        equal: comparison_unit,
        logical_and: comparison_unit,
        logical_or: comparison_unit,
        logical_xor: comparison_unit,
        logical_not: return_without_unit,
        maximum: preserve_units,
        minimum: preserve_units,
        fmax: preserve_units,
        fmin: preserve_units,
        isreal: return_without_unit,
        iscomplex: return_without_unit,
        isfinite: return_without_unit,
        isinf: return_without_unit,
        isnan: return_without_unit,
        signbit: return_without_unit,
        copysign: passthrough_unit,
        nextafter: preserve_units,
        modf: passthrough_unit,
        ldexp: bitop_units,
        frexp: return_without_unit,
        floor: passthrough_unit,
        ceil: passthrough_unit,
        trunc: passthrough_unit,
        spacing: passthrough_unit,
    }

    __array_priority__ = 2.0

    def __new__(cls, input_array, input_units=None, registry=None, dtype=None,
                bypass_validation=False):
        if dtype is None:
            dtype = getattr(input_array, 'dtype', np.float64)
        if bypass_validation is True:
            obj = np.asarray(input_array, dtype=dtype).view(cls)
            obj.units = input_units
            if registry is not None:
                obj.units.registry = registry
            return obj
        if input_array is NotImplemented:
            return input_array.view(cls)
        if registry is None and isinstance(input_units, (str, bytes)):
            if input_units.startswith('code_'):
                raise UnitParseError(
                    "Code units used without referring to a dataset. \n"
                    "Perhaps you meant to do something like this instead: \n"
                    "ds.arr(%s, \"%s\")" % (input_array, input_units)
                    )
        if isinstance(input_array, YTArray):
            ret = input_array.view(cls)
            if input_units is None:
                if registry is None:
                    pass
                else:
                    units = Unit(str(input_array.units), registry=registry)
                    ret.units = units
            elif isinstance(input_units, Unit):
                ret.units = input_units
            else:
                ret.units = Unit(input_units, registry=registry)
            return ret
        elif isinstance(input_array, np.ndarray):
            pass
        elif iterable(input_array) and input_array:
            if isinstance(input_array[0], YTArray):
                return YTArray(np.array(input_array, dtype=dtype),
                               input_array[0].units, registry=registry)

        # Input array is an already formed ndarray instance
        # We first cast to be our class type

        obj = np.asarray(input_array, dtype=dtype).view(cls)

        # Check units type
        if input_units is None:
            # Nothing provided. Make dimensionless...
            units = Unit()
        elif isinstance(input_units, Unit):
            if registry and registry is not input_units.registry:
                units = Unit(str(input_units), registry=registry)
            else:
                units = input_units
        else:
            # units kwarg set, but it's not a Unit object.
            # don't handle all the cases here, let the Unit class handle if
            # it's a str.
            units = Unit(input_units, registry=registry)

        # Attach the units
        obj.units = units

        return obj

    def __array_finalize__(self, obj):
        """

        """
        if obj is None and hasattr(self, 'units'):
            return
        self.units = getattr(obj, 'units', NULL_UNIT)

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

    def convert_to_units(self, units):
        """
        Convert the array and units to the given units.

        Parameters
        ----------
        units : Unit object or str
            The units you want to convert to.

        """
        new_units = _unit_repr_check_same(self.units, units)
        (conversion_factor, offset) = self.units.get_conversion_factor(new_units)

        self.units = new_units
        values = self.d
        values *= conversion_factor

        if offset:
            np.subtract(self, offset*self.uq, self)

        return self

    def convert_to_base(self, unit_system="cgs"):
        """
        Convert the array and units to the equivalent base units in
        the specified unit system.

        Parameters
        ----------
        unit_system : string, optional
            The unit system to be used in the conversion. If not specified,
            the default base units of cgs are used.

        Examples
        --------
        >>> E = YTQuantity(2.5, "erg/s")
        >>> E.convert_to_base(unit_system="galactic")
        """
        return self.convert_to_units(self.units.get_base_equivalent(unit_system))

    def convert_to_cgs(self):
        """
        Convert the array and units to the equivalent cgs units.

        """
        return self.convert_to_units(self.units.get_cgs_equivalent())

    def convert_to_mks(self):
        """
        Convert the array and units to the equivalent mks units.

        """
        return self.convert_to_units(self.units.get_mks_equivalent())

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
        new_units = _unit_repr_check_same(self.units, units)
        (conversion_factor, offset) = self.units.get_conversion_factor(new_units)

        new_array = type(self)(self.ndview * conversion_factor, new_units)

        if offset:
            np.subtract(new_array, offset*new_array.uq, new_array)

        return new_array

    def to(self, units):
        """
        An alias for YTArray.in_units().

        See the docstrings of that function for details.
        """
        return self.in_units(units)

    def in_base(self, unit_system="cgs"):
        """
        Creates a copy of this array with the data in the specified unit system,
        and returns it in that system's base units.

        Parameters
        ----------
        unit_system : string, optional
            The unit system to be used in the conversion. If not specified,
            the default base units of cgs are used.

        Examples
        --------
        >>> E = YTQuantity(2.5, "erg/s")
        >>> E_new = E.in_base(unit_system="galactic")
        """
        return self.in_units(self.units.get_base_equivalent(unit_system))

    def in_cgs(self):
        """
        Creates a copy of this array with the data in the equivalent cgs units,
        and returns it.

        Returns
        -------
        Quantity object with data converted to cgs units.

        """
        return self.in_units(self.units.get_cgs_equivalent())

    def in_mks(self):
        """
        Creates a copy of this array with the data in the equivalent mks units,
        and returns it.

        Returns
        -------
        Quantity object with data converted to mks units.

        """
        return self.in_units(self.units.get_mks_equivalent())

    def to_equivalent(self, unit, equiv, **kwargs):
        """
        Convert a YTArray or YTQuantity to an equivalent, e.g., something that is
        related by only a constant factor but not in the same units.

        Parameters
        ----------
        unit : string
            The unit that you wish to convert to.
        equiv : string
            The equivalence you wish to use. To see which equivalencies are
            supported for this unitful quantity, try the
            :meth:`list_equivalencies` method.

        Examples
        --------
        >>> a = yt.YTArray(1.0e7,"K")
        >>> a.to_equivalent("keV", "thermal")
        """
        conv_unit = Unit(unit, registry=self.units.registry)
        this_equiv = equivalence_registry[equiv]()
        oneway_or_equivalent = (
            conv_unit.has_equivalent(equiv) or this_equiv._one_way)
        if self.has_equivalent(equiv) and oneway_or_equivalent:
            new_arr = this_equiv.convert(
                self, conv_unit.dimensions, **kwargs)
            if isinstance(new_arr, tuple):
                try:
                    return YTArray(new_arr[0], new_arr[1]).in_units(unit)
                except YTUnitConversionError:
                    raise YTInvalidUnitEquivalence(equiv, self.units, unit)
            else:
                return new_arr.in_units(unit)
        else:
            raise YTInvalidUnitEquivalence(equiv, self.units, unit)

    def list_equivalencies(self):
        """
        Lists the possible equivalencies associated with this YTArray or
        YTQuantity.
        """
        self.units.list_equivalencies()

    def has_equivalent(self, equiv):
        """
        Check to see if this YTArray or YTQuantity has an equivalent unit in
        *equiv*.
        """
        return self.units.has_equivalent(equiv)

    def ndarray_view(self):
        """
        Returns a view into the array, but as an ndarray rather than ytarray.

        Returns
        -------
        View of this array's data.
        """
        return self.view(np.ndarray)

    def to_ndarray(self):
        """
        Creates a copy of this array with the unit information stripped

        """
        return np.array(self)

    @classmethod
    def from_astropy(cls, arr, unit_registry=None):
        """
        Convert an AstroPy "Quantity" to a YTArray or YTQuantity.

        Parameters
        ----------
        arr : AstroPy Quantity
            The Quantity to convert from.
        unit_registry : yt UnitRegistry, optional
            A yt unit registry to use in the conversion. If one is not
            supplied, the default one will be used.
        """
        # Converting from AstroPy Quantity
        u = arr.unit
        ap_units = []
        for base, exponent in zip(u.bases, u.powers):
            unit_str = base.to_string()
            # we have to do this because AstroPy is silly and defines
            # hour as "h"
            if unit_str == "h": unit_str = "hr"
            ap_units.append("%s**(%s)" % (unit_str, Rational(exponent)))
        ap_units = "*".join(ap_units)
        if isinstance(arr.value, np.ndarray):
            return YTArray(arr.value, ap_units, registry=unit_registry)
        else:
            return YTQuantity(arr.value, ap_units, registry=unit_registry)


    def to_astropy(self, **kwargs):
        """
        Creates a new AstroPy quantity with the same unit information.
        """
        if _astropy.units is None:
            raise ImportError("You don't have AstroPy installed, so you can't convert to " +
                              "an AstroPy quantity.")
        return self.value*_astropy.units.Unit(str(self.units), **kwargs)

    @classmethod
    def from_pint(cls, arr, unit_registry=None):
        """
        Convert a Pint "Quantity" to a YTArray or YTQuantity.

        Parameters
        ----------
        arr : Pint Quantity
            The Quantity to convert from.
        unit_registry : yt UnitRegistry, optional
            A yt unit registry to use in the conversion. If one is not
            supplied, the default one will be used.

        Examples
        --------
        >>> from pint import UnitRegistry
        >>> import numpy as np
        >>> ureg = UnitRegistry()
        >>> a = np.random.random(10)
        >>> b = ureg.Quantity(a, "erg/cm**3")
        >>> c = yt.YTArray.from_pint(b)
        """
        p_units = []
        for base, exponent in arr._units.items():
            bs = convert_pint_units(base)
            p_units.append("%s**(%s)" % (bs, Rational(exponent)))
        p_units = "*".join(p_units)
        if isinstance(arr.magnitude, np.ndarray):
            return YTArray(arr.magnitude, p_units, registry=unit_registry)
        else:
            return YTQuantity(arr.magnitude, p_units, registry=unit_registry)

    def to_pint(self, unit_registry=None):
        """
        Convert a YTArray or YTQuantity to a Pint Quantity.

        Parameters
        ----------
        arr : YTArray or YTQuantity
            The unitful quantity to convert from.
        unit_registry : Pint UnitRegistry, optional
            The Pint UnitRegistry to use in the conversion. If one is not
            supplied, the default one will be used. NOTE: This is not
            the same as a yt UnitRegistry object.

        Examples
        --------
        >>> a = YTQuantity(4.0, "cm**2/s")
        >>> b = a.to_pint()
        """
        from pint import UnitRegistry
        if unit_registry is None:
            unit_registry = UnitRegistry()
        powers_dict = self.units.expr.as_powers_dict()
        units = []
        for unit, pow in powers_dict.items():
            # we have to do this because Pint doesn't recognize
            # "yr" as "year"
            if str(unit).endswith("yr") and len(str(unit)) in [2,3]:
                unit = str(unit).replace("yr","year")
            units.append("%s**(%s)" % (unit, Rational(pow)))
        units = "*".join(units)
        return unit_registry.Quantity(self.value, units)

    #
    # End unit conversion methods
    #

    def write_hdf5(self, filename, dataset_name=None, info=None, group_name=None):
        r"""Writes a YTArray to hdf5 file.

        Parameters
        ----------
        filename: string
            The filename to create and write a dataset to

        dataset_name: string
            The name of the dataset to create in the file.

        info: dictionary
            A dictionary of supplementary info to write to append as attributes
            to the dataset.
            
        group_name: string
            An optional group to write the arrays to. If not specified, the arrays
            are datasets at the top level by default.

        Examples
        --------
        >>> a = YTArray([1,2,3], 'cm')
        >>> myinfo = {'field':'dinosaurs', 'type':'field_data'}
        >>> a.write_hdf5('test_array_data.h5', dataset_name='dinosaurs',
        ...              info=myinfo)
        """
        import h5py
        from yt.extern.six.moves import cPickle as pickle
        if info is None:
            info = {}

        info['units'] = str(self.units)
        info['unit_registry'] = np.void(pickle.dumps(self.units.registry.lut))

        if dataset_name is None:
            dataset_name = 'array_data'

        f = h5py.File(filename)
        if group_name is not None:
            if group_name in f:
                g = f[group_name]
            else:
                g = f.create_group(group_name)
        else:
            g = f
        if dataset_name in g.keys():
            d = g[dataset_name]
            # Overwrite without deleting if we can get away with it.
            if d.shape == self.shape and d.dtype == self.dtype:
                d[:] = self
                for k in d.attrs.keys():
                    del d.attrs[k]
            else:
                del f[dataset_name]
                d = g.create_dataset(dataset_name, data=self)
        else:
            d = g.create_dataset(dataset_name, data=self)

        for k, v in info.items():
            d.attrs[k] = v
        f.close()

    @classmethod
    def from_hdf5(cls, filename, dataset_name=None, group_name=None):
        r"""Attempts read in and convert a dataset in an hdf5 file into a YTArray.

        Parameters
        ----------
        filename: string
        The filename to of the hdf5 file.

        dataset_name: string
            The name of the dataset to read from.  If the dataset has a units
            attribute, attempt to infer units as well.

        group_name: string
            An optional group to read the arrays from. If not specified, the arrays
            are datasets at the top level by default.

        """
        import h5py
        from yt.extern.six.moves import cPickle as pickle

        if dataset_name is None:
            dataset_name = 'array_data'

        f = h5py.File(filename)
        if group_name is not None:
            g = f[group_name]
        else:
            g = f
        dataset = g[dataset_name]
        data = dataset[:]
        units = dataset.attrs.get('units', '')
        if 'unit_registry' in dataset.attrs.keys():
            unit_lut = pickle.loads(dataset.attrs['unit_registry'].tostring())
        else:
            unit_lut = None
        f.close()
        registry = UnitRegistry(lut=unit_lut, add_default_symbols=False)
        return cls(data, units, registry=registry)

    #
    # Start convenience methods
    #

    @property
    def value(self):
        """Get a copy of the array data as a numpy ndarray"""
        return np.array(self)

    v = value

    @property
    def ndview(self):
        """Get a view of the array data."""
        return self.ndarray_view()

    d = ndview

    @property
    def unit_quantity(self):
        """Get a YTQuantity with the same unit as this array and a value of 1.0"""
        return YTQuantity(1.0, self.units)

    uq = unit_quantity

    @property
    def unit_array(self):
        """Get a YTArray filled with ones with the same unit and shape as this array"""
        return np.ones_like(self)

    ua = unit_array

    #
    # Start operation methods
    #

    def __add__(self, right_object):
        """
        Add this ytarray to the object on the right of the `+` operator. Must
        check for the correct (same dimension) units.

        """
        ro = sanitize_units_add(self, right_object, "addition")
        return super(YTArray, self).__add__(ro)

    def __radd__(self, left_object):
        """ See __add__. """
        lo = sanitize_units_add(self, left_object, "addition")
        return super(YTArray, self).__radd__(lo)

    def __iadd__(self, other):
        """ See __add__. """
        oth = sanitize_units_add(self, other, "addition")
        np.add(self, oth, out=self)
        return self

    def __sub__(self, right_object):
        """
        Subtract the object on the right of the `-` from this ytarray. Must
        check for the correct (same dimension) units.

        """
        ro = sanitize_units_add(self, right_object, "subtraction")
        return super(YTArray, self).__sub__(ro)

    def __rsub__(self, left_object):
        """ See __sub__. """
        lo = sanitize_units_add(self, left_object, "subtraction")
        return super(YTArray, self).__rsub__(lo)

    def __isub__(self, other):
        """ See __sub__. """
        oth = sanitize_units_add(self, other, "subtraction")
        np.subtract(self, oth, out=self)
        return self

    def __neg__(self):
        """ Negate the data. """
        return super(YTArray, self).__neg__()

    def __pos__(self):
        """ Posify the data. """
        return type(self)(super(YTArray, self).__pos__(), self.units)

    def __mul__(self, right_object):
        """
        Multiply this YTArray by the object on the right of the `*` operator.
        The unit objects handle being multiplied.

        """
        ro = sanitize_units_mul(self, right_object)
        return super(YTArray, self).__mul__(ro)

    def __rmul__(self, left_object):
        """ See __mul__. """
        lo = sanitize_units_mul(self, left_object)
        return super(YTArray, self).__rmul__(lo)

    def __imul__(self, other):
        """ See __mul__. """
        oth = sanitize_units_mul(self, other)
        np.multiply(self, oth, out=self)
        return self

    def __div__(self, right_object):
        """
        Divide this YTArray by the object on the right of the `/` operator.

        """
        ro = sanitize_units_mul(self, right_object)
        return super(YTArray, self).__div__(ro)

    def __rdiv__(self, left_object):
        """ See __div__. """
        lo = sanitize_units_mul(self, left_object)
        return super(YTArray, self).__rdiv__(lo)

    def __idiv__(self, other):
        """ See __div__. """
        oth = sanitize_units_mul(self, other)
        np.divide(self, oth, out=self)
        return self

    def __truediv__(self, right_object):
        ro = sanitize_units_mul(self, right_object)
        return super(YTArray, self).__truediv__(ro)

    def __rtruediv__(self, left_object):
        """ See __div__. """
        lo = sanitize_units_mul(self, left_object)
        return super(YTArray, self).__rtruediv__(lo)

    def __itruediv__(self, other):
        """ See __div__. """
        oth = sanitize_units_mul(self, other)
        np.true_divide(self, oth, out=self)
        return self

    def __floordiv__(self, right_object):
        ro = sanitize_units_mul(self, right_object)
        return super(YTArray, self).__floordiv__(ro)

    def __rfloordiv__(self, left_object):
        """ See __div__. """
        lo = sanitize_units_mul(self, left_object)
        return super(YTArray, self).__rfloordiv__(lo)

    def __ifloordiv__(self, other):
        """ See __div__. """
        oth = sanitize_units_mul(self, other)
        np.floor_divide(self, oth, out=self)
        return self

    def __or__(self, right_object):
        return super(YTArray, self).__or__(right_object)

    def __ror__(self, left_object):
        return super(YTArray, self).__ror__(left_object)

    def __ior__(self, other):
        np.bitwise_or(self, other, out=self)
        return self

    def __xor__(self, right_object):
        return super(YTArray, self).__xor__(right_object)

    def __rxor__(self, left_object):
        return super(YTArray, self).__rxor__(left_object)

    def __ixor__(self, other):
        np.bitwise_xor(self, other, out=self)
        return self

    def __and__(self, right_object):
        return super(YTArray, self).__and__(right_object)

    def __rand__(self, left_object):
        return super(YTArray, self).__rand__(left_object)

    def __iand__(self, other):
        np.bitwise_and(self, other, out=self)
        return self

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
                raise YTUnitOperationError('power', power.unit)

        # Work around a sympy issue (I think?)
        #
        # If I don't do this, super(YTArray, self).__pow__ returns a YTArray
        # with a unit attribute set to the sympy expression 1/1 rather than a
        # dimensionless Unit object.
        if self.units.is_dimensionless and power == -1:
            ret = super(YTArray, self).__pow__(power)
            return type(self)(ret, input_units='')

        return super(YTArray, self).__pow__(power)

    def __abs__(self):
        """ Return a YTArray with the abs of the data. """
        return super(YTArray, self).__abs__()

    def sqrt(self):
        """
        Return sqrt of this YTArray. We take the sqrt for the array and use
        take the 1/2 power of the units.

        """
        return type(self)(super(YTArray, self).sqrt(),
                          input_units=self.units**0.5)

    #
    # Start comparison operators.
    #

    def __lt__(self, other):
        """ Test if this is less than the object on the right. """
        # converts if possible
        oth = validate_comparison_units(self, other, 'less_than')
        return super(YTArray, self).__lt__(oth)

    def __le__(self, other):
        """ Test if this is less than or equal to the object on the right. """
        oth = validate_comparison_units(self, other, 'less_than or equal')
        return super(YTArray, self).__le__(oth)

    def __eq__(self, other):
        """ Test if this is equal to the object on the right. """
        # Check that other is a YTArray.
        if other is None:
            # self is a YTArray, so it can't be None.
            return False
        oth = validate_comparison_units(self, other, 'equal')
        return super(YTArray, self).__eq__(oth)

    def __ne__(self, other):
        """ Test if this is not equal to the object on the right. """
        # Check that the other is a YTArray.
        if other is None:
            return True
        oth = validate_comparison_units(self, other, 'not equal')
        return super(YTArray, self).__ne__(oth)

    def __ge__(self, other):
        """ Test if this is greater than or equal to other. """
        # Check that the other is a YTArray.
        oth = validate_comparison_units(self, other, 'greater than or equal')
        return super(YTArray, self).__ge__(oth)

    def __gt__(self, other):
        """ Test if this is greater than the object on the right. """
        # Check that the other is a YTArray.
        oth = validate_comparison_units(self, other, 'greater than')
        return super(YTArray, self).__gt__(oth)

    #
    # End comparison operators
    #

    #
    # Begin reduction operators
    #

    @return_arr
    def prod(self, axis=None, dtype=None, out=None):
        if axis is not None:
            units = self.units**self.shape[axis]
        else:
            units = self.units**self.size
        return super(YTArray, self).prod(axis, dtype, out), units

    @return_arr
    def mean(self, axis=None, dtype=None, out=None):
        return super(YTArray, self).mean(axis, dtype, out), self.units

    @return_arr
    def sum(self, axis=None, dtype=None, out=None):
        return super(YTArray, self).sum(axis, dtype, out), self.units

    @return_arr
    def dot(self, b, out=None):
        return super(YTArray, self).dot(b), self.units*b.units

    @return_arr
    def std(self, axis=None, dtype=None, out=None, ddof=0):
        return super(YTArray, self).std(axis, dtype, out, ddof), self.units

    def __getitem__(self, item):
        ret = super(YTArray, self).__getitem__(item)
        if ret.shape == ():
            return YTQuantity(ret, self.units, bypass_validation=True)
        else:
            return ret

    def __array_wrap__(self, out_arr, context=None):
        ret = super(YTArray, self).__array_wrap__(out_arr, context)
        if isinstance(ret, YTQuantity) and ret.shape != ():
            ret = ret.view(YTArray)
        if context is None:
            if ret.shape == ():
                return ret[()]
            else:
                return ret
        elif context[0] in unary_operators:
            u = getattr(context[1][0], 'units', None)
            if u is None:
                u = NULL_UNIT
            if u.dimensions is angle and context[0] in trigonometric_operators:
                out_arr = context[0](
                    context[1][0].in_units('radian').view(np.ndarray))
            unit = self._ufunc_registry[context[0]](u)
            ret_class = type(self)
        elif context[0] in binary_operators:
            oper1 = coerce_iterable_units(context[1][0])
            oper2 = coerce_iterable_units(context[1][1])
            cls1 = type(oper1)
            cls2 = type(oper2)
            unit1 = getattr(oper1, 'units', None)
            unit2 = getattr(oper2, 'units', None)
            ret_class = get_binary_op_return_class(cls1, cls2)
            if unit1 is None:
                unit1 = Unit(registry=getattr(unit2, 'registry', None))
            if unit2 is None and context[0] is not power:
                unit2 = Unit(registry=getattr(unit1, 'registry', None))
            elif context[0] is power:
                unit2 = oper2
                if isinstance(unit2, np.ndarray):
                    if isinstance(unit2, YTArray):
                        if unit2.units.is_dimensionless:
                            pass
                        else:
                            raise YTUnitOperationError(context[0], unit1, unit2)
                    unit2 = 1.0
            unit_operator = self._ufunc_registry[context[0]]
            if unit_operator in (preserve_units, comparison_unit, arctan2_unit):
                # Allow comparisons, addition, and subtraction with
                # dimensionless quantities or arrays filled with zeros.
                u1d = unit1.is_dimensionless
                u2d = unit2.is_dimensionless
                if unit1 != unit2:
                    any_nonzero = [np.any(oper1), np.any(oper2)]
                    if any_nonzero[0] is np.bool_(False):
                        unit1 = unit2
                    elif any_nonzero[1] is np.bool_(False):
                        unit2 = unit1
                    elif not any([u1d, u2d]):
                        if not unit1.same_dimensions_as(unit2):
                            raise YTUnitOperationError(context[0], unit1, unit2)
                        else:
                            raise YTUfuncUnitError(context[0], unit1, unit2)
            unit = unit_operator(unit1, unit2)
            if unit_operator in (multiply_units, divide_units):
                if unit.is_dimensionless and unit.base_value != 1.0:
                    if not unit1.is_dimensionless:
                        if unit1.dimensions == unit2.dimensions:
                            np.multiply(out_arr.view(np.ndarray),
                                        unit.base_value, out=out_arr)
                            unit = Unit(registry=unit.registry)
        else:
            raise RuntimeError("Support for the %s ufunc has not been added "
                               "to YTArray." % str(context[0]))
        if unit is None:
            out_arr = np.array(out_arr, copy=False)
            return out_arr
        out_arr.units = unit
        if out_arr.size == 1:
            return YTQuantity(np.array(out_arr), unit)
        else:
            if ret_class is YTQuantity:
                # This happens if you do ndarray * YTQuantity. Explicitly
                # casting to YTArray avoids creating a YTQuantity with size > 1
                return YTArray(np.array(out_arr), unit)
            return ret_class(np.array(out_arr, copy=False), unit)


    def __reduce__(self):
        """Pickle reduction method

        See the documentation for the standard library pickle module:
        http://docs.python.org/2/library/pickle.html

        Unit metadata is encoded in the zeroth element of third element of the
        returned tuple, itself a tuple used to restore the state of the ndarray.
        This is always defined for numpy arrays.
        """
        np_ret = super(YTArray, self).__reduce__()
        obj_state = np_ret[2]
        unit_state = (((str(self.units), self.units.registry.lut),) + obj_state[:],)
        new_ret = np_ret[:2] + unit_state + np_ret[3:]
        return new_ret

    def __setstate__(self, state):
        """Pickle setstate method

        This is called inside pickle.read() and restores the unit data from the
        metadata extracted in __reduce__ and then serialized by pickle.
        """
        super(YTArray, self).__setstate__(state[1:])
        unit, lut = state[0]
        # need to fix up the lut if the pickle was saved prior to PR #1728
        # when the pickle format changed
        if len(lut['m']) == 2:
            lut.update(default_unit_symbol_lut)
            for k, v in [(k, v) for k, v in lut.items() if len(v) == 2]:
                lut[k] = v + (0.0, r'\rm{' + k.replace('_', '\ ') + '}')
        registry = UnitRegistry(lut=lut, add_default_symbols=False)
        self.units = Unit(unit, registry=registry)

    def __deepcopy__(self, memodict=None):
        """copy.deepcopy implementation

        This is necessary for stdlib deepcopy of arrays and quantities.
        """
        if memodict is None:
            memodict = {}
        ret = super(YTArray, self).__deepcopy__(memodict)
        return type(self)(ret, copy.deepcopy(self.units))

class YTQuantity(YTArray):
    """
    A scalar associated with a unit.

    Parameters
    ----------

    input_scalar : an integer or floating point scalar
        The scalar to attach units to
    input_units : String unit specification, unit symbol object, or astropy units
        The units of the quantity. Powers must be specified using python syntax
        (cm**3, not cm^3).
    registry : A UnitRegistry object
        The registry to create units from. If input_units is already associated
        with a unit registry and this is specified, this will be used instead of
        the registry associated with the unit object.
    dtype : string or NumPy dtype object
        The dtype of the array data.

    Examples
    --------

    >>> from yt import YTQuantity
    >>> a = YTQuantity(1, 'cm')
    >>> b = YTQuantity(2, 'm')
    >>> a + b
    201.0 cm
    >>> b + a
    2.01 m

    NumPy ufuncs will pass through units where appropriate.

    >>> import numpy as np
    >>> a = YTQuantity(12, 'g/cm**3')
    >>> np.ones_like(a)
    1 g/cm**3

    and strip them when it would be annoying to deal with them.

    >>> print(np.log10(a))
    1.07918124605

    YTQuantity is tightly integrated with yt datasets:

    >>> import yt
    >>> ds = yt.load('IsolatedGalaxy/galaxy0030/galaxy0030')
    >>> a = ds.quan(5, 'code_length')
    >>> a.in_cgs()
    1.543e+25 cm

    This is equivalent to:

    >>> b = YTQuantity(5, 'code_length', registry=ds.unit_registry)
    >>> np.all(a == b)
    True

    """
    def __new__(cls, input_scalar, input_units=None, registry=None,
                dtype=np.float64, bypass_validation=False):
        if not isinstance(input_scalar, (numeric_type, np.number, np.ndarray)):
            raise RuntimeError("YTQuantity values must be numeric")
        ret = YTArray.__new__(cls, input_scalar, input_units, registry,
                              dtype=dtype, bypass_validation=bypass_validation)
        if ret.size > 1:
            raise RuntimeError("YTQuantity instances must be scalars")
        return ret

    def __repr__(self):
        return str(self)

def validate_numpy_wrapper_units(v, arrs):
    if not any(isinstance(a, YTArray) for a in arrs):
        return v
    if not all(isinstance(a, YTArray) for a in arrs):
        raise RuntimeError("Not all of your arrays are YTArrays.")
    a1 = arrs[0]
    if not all(a.units == a1.units for a in arrs[1:]):
        raise RuntimeError("Your arrays must have identical units.")
    v.units = a1.units
    return v

def uconcatenate(arrs, axis=0):
    """Concatenate a sequence of arrays.

    This wrapper around numpy.concatenate preserves units. All input arrays must
    have the same units.  See the documentation of numpy.concatenate for full
    details.

    Examples
    --------
    >>> A = yt.YTArray([1, 2, 3], 'cm')
    >>> B = yt.YTArray([2, 3, 4], 'cm')
    >>> uconcatenate((A, B))
    YTArray([ 1., 2., 3., 2., 3., 4.]) cm

    """
    v = np.concatenate(arrs, axis=axis)
    v = validate_numpy_wrapper_units(v, arrs)
    return v

def ucross(arr1,arr2, registry=None):
    """Applies the cross product to two YT arrays.

    This wrapper around numpy.cross preserves units.
    See the documentation of numpy.cross for full
    details.
    """

    v = np.cross(arr1,arr2)
    units = arr1.units * arr2.units
    arr = YTArray(v,units, registry=registry)
    return arr

def uintersect1d(arr1, arr2, assume_unique=False):
    """Find the sorted unique elements of the two input arrays.

    A wrapper around numpy.intersect1d that preserves units.  All input arrays
    must have the same units.  See the documentation of numpy.intersect1d for
    full details.

    Examples
    --------
    >>> A = yt.YTArray([1, 2, 3], 'cm')
    >>> B = yt.YTArray([2, 3, 4], 'cm')
    >>> uintersect1d(A, B)
    YTArray([ 2., 3.]) cm

    """
    v = np.intersect1d(arr1, arr2, assume_unique=assume_unique)
    v = validate_numpy_wrapper_units(v, [arr1, arr2])
    return v

def uunion1d(arr1, arr2):
    """Find the union of two arrays.

    A wrapper around numpy.intersect1d that preserves units.  All input arrays
    must have the same units.  See the documentation of numpy.intersect1d for
    full details.

    Examples
    --------
    >>> A = yt.YTArray([1, 2, 3], 'cm')
    >>> B = yt.YTArray([2, 3, 4], 'cm')
    >>> uunion1d(A, B)
    YTArray([ 1., 2., 3., 4.]) cm

    """
    v = np.union1d(arr1, arr2)
    v = validate_numpy_wrapper_units(v, [arr1, arr2])
    return v

def unorm(data):
    """Matrix or vector norm that preserves units

    This is a wrapper around np.linalg.norm that preserves units.
    """
    return YTArray(np.linalg.norm(data), data.units)

def uvstack(arrs):
    """Stack arrays in sequence vertically (row wise) while preserving units

    This is a wrapper around np.vstack that preserves units.
    """
    v = np.vstack(arrs)
    v = validate_numpy_wrapper_units(v, arrs)
    return v

def uhstack(arrs):
    """Stack arrays in sequence horizontally (column wise) while preserving units

    This is a wrapper around np.hstack that preserves units.
    """
    v = np.hstack(arrs)
    v = validate_numpy_wrapper_units(v, arrs)
    return v

def array_like_field(data, x, field):
    field = data._determine_fields(field)[0]
    if isinstance(field, tuple):
        units = data.ds._get_field_info(field[0],field[1]).output_units
    else:
        units = data.ds._get_field_info(field).output_units
    if isinstance(x, YTArray):
        arr = copy.deepcopy(x)
        arr.convert_to_units(units)
        return arr
    if isinstance(x, np.ndarray):
        return data.ds.arr(x, units)
    else:
        return data.ds.quan(x, units)

def get_binary_op_return_class(cls1, cls2):
    if cls1 is cls2:
        return cls1
    if cls1 is np.ndarray or issubclass(cls1, (numeric_type, np.number, list, tuple)):
        return cls2
    if cls2 is np.ndarray or issubclass(cls2, (numeric_type, np.number, list, tuple)):
        return cls1
    if issubclass(cls1, YTQuantity):
        return cls2
    if issubclass(cls2, YTQuantity):
        return cls1
    if issubclass(cls1, cls2):
        return cls1
    if issubclass(cls2, cls1):
        return cls2
    else:
        raise RuntimeError("Undefined operation for a YTArray subclass. "
                           "Received operand types (%s) and (%s)" % (cls1, cls2))

def loadtxt(fname, dtype='float', delimiter='\t', usecols=None, comments='#'):
    r"""
    Load YTArrays with unit information from a text file. Each row in the
    text file must have the same number of values.

    Parameters
    ----------
    fname : str
        Filename to read.
    dtype : data-type, optional
        Data-type of the resulting array; default: float.
    delimiter : str, optional
        The string used to separate values.  By default, this is any
        whitespace.
    usecols : sequence, optional
        Which columns to read, with 0 being the first.  For example,
        ``usecols = (1,4,5)`` will extract the 2nd, 5th and 6th columns.
        The default, None, results in all columns being read.
    comments : str, optional
        The character used to indicate the start of a comment;
        default: '#'.

    Examples
    --------
    >>> temp, velx = yt.loadtxt("sphere.dat", usecols=(1,2), delimiter="\t")
    """
    f = open(fname, 'r')
    next_one = False
    units = []
    num_cols = -1
    for line in f.readlines():
        words = line.strip().split()
        if len(words) == 0:
            continue
        if line[0] == comments:
            if next_one:
                units = words[1:]
            if len(words) == 2 and words[1] == "Units":
                next_one = True
        else:
            # Here we catch the first line of numbers
            try:
                col_words = line.strip().split(delimiter)
                for word in col_words:
                    float(word)
                num_cols = len(col_words)
                break
            except ValueError:
                mylog.warning("Unrecognized character at beginning of line: \"%s\"." % line[0])
    f.close()
    if len(units) != num_cols:
        mylog.warning("Malformed or incomplete units header. Arrays will be "
                      "dimensionless!")
        units = ["dimensionless"]*num_cols
    arrays = np.loadtxt(fname, dtype=dtype, comments=comments,
                        delimiter=delimiter, converters=None,
                        unpack=True, usecols=usecols, ndmin=0)
    if usecols is not None:
        units = [units[col] for col in usecols]
    mylog.info("Array units: %s" % ", ".join(units))
    return tuple([YTArray(arr, unit) for arr, unit in zip(arrays, units)])

def savetxt(fname, arrays, fmt='%.18e', delimiter='\t', header='',
            footer='', comments='#'):
    r"""
    Write YTArrays with unit information to a text file.

    Parameters
    ----------
    fname : str
        The file to write the YTArrays to.
    arrays : list of YTArrays or single YTArray
        The array(s) to write to the file.
    fmt : str or sequence of strs, optional
        A single format (%10.5f), or a sequence of formats.
    delimiter : str, optional
        String or character separating columns.
    header : str, optional
        String that will be written at the beginning of the file, before the
        unit header.
    footer : str, optional
        String that will be written at the end of the file.
    comments : str, optional
        String that will be prepended to the ``header`` and ``footer`` strings,
        to mark them as comments. Default: '# ', as expected by e.g.
        ``yt.loadtxt``.

    Examples
    --------
    >>> sp = ds.sphere("c", (100,"kpc"))
    >>> a = sp["density"]
    >>> b = sp["temperature"]
    >>> c = sp["velocity_x"]
    >>> yt.savetxt("sphere.dat", [a,b,c], header='My sphere stuff', delimiter="\t")
    """
    if not isinstance(arrays, list):
        arrays = [arrays]
    units = []
    for array in arrays:
        if hasattr(array, "units"):
            units.append(str(array.units))
        else:
            units.append("dimensionless")
    if header != '':
        header += '\n'
    header += " Units\n " + '\t'.join(units)
    np.savetxt(fname, np.transpose(arrays), header=header,
               fmt=fmt, delimiter=delimiter, footer=footer,
               newline='\n', comments=comments)
