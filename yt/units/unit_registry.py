"""
A registry for units that can be added to and modified.


"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import json
import re

from distutils.version import LooseVersion
from yt.units.unit_lookup_table import \
    default_unit_symbol_lut
from yt.utilities.lib.fnv_hash import fnv_hash
from yt.extern import six
from sympy import \
    sympify, \
    srepr, \
    __version__ as sympy_version

SYMPY_VERSION = LooseVersion(sympy_version)

def positive_symbol_replacer(match):
    return match.group().replace(')\')', ')\', positive=True)')

class SymbolNotFoundError(Exception):
    pass

class UnitParseError(Exception):
    pass

class UnitRegistry:
    """A registry for unit symbols"""
    def __init__(self, add_default_symbols=True, lut=None):
        if lut:
            self.lut = lut
        else:
            self.lut = {}
        self.unit_objs = {}

        if add_default_symbols:
            self.lut.update(default_unit_symbol_lut)

    def __getitem__(self, key):
        return self.lut[key]

    def __contains__(self, item):
        return item in self.lut

    _unit_system_id = None
    @property
    def unit_system_id(self):
        """
        This is a unique identifier for the unit registry created 
        from a FNV hash. It is needed to register a dataset's code
        unit system in the unit system registry.
        """
        if self._unit_system_id is None:
            hash_data = bytearray()
            for k, v in self.lut.items():
                hash_data.extend(k.encode('ascii'))
                hash_data.extend(repr(v).encode('ascii'))
            self._unit_system_id = "us_%d" % fnv_hash(hash_data)
        return self._unit_system_id

    def add(self, symbol, base_value, dimensions, tex_repr=None, offset=None):
        """
        Add a symbol to this registry.

        """
        from yt.units.unit_object import validate_dimensions

        # Validate
        if not isinstance(base_value, float):
            raise UnitParseError("base_value (%s) must be a float, got a %s."
                                 % (base_value, type(base_value)))

        if offset is not None:
            if not isinstance(offset, float):
                raise UnitParseError(
                    "offset value (%s) must be a float, got a %s."
                    % (offset, type(offset)))
        else:
            offset = 0.0

        validate_dimensions(dimensions)

        if tex_repr is None:
            # make educated guess that will look nice in most cases
            tex_repr = r"\rm{" + symbol.replace('_', '\ ') + "}"

        # Add to lut
        self.lut.update({symbol: (base_value, dimensions, offset, tex_repr)})

    def remove(self, symbol):
        """
        Remove the entry for the unit matching `symbol`.

        """
        if symbol not in self.lut:
            raise SymbolNotFoundError(
                "Tried to remove the symbol '%s', but it does not exist" \
                "in this registry." % symbol)

        del self.lut[symbol]
        if symbol in self.unit_objs:
            del self.unit_objs[symbol]

    def modify(self, symbol, base_value):
        """
        Change the base value of a unit symbol.  Useful for adjusting code units
        after parsing parameters.

        """
        if symbol not in self.lut:
            raise SymbolNotFoundError(
                "Tried to modify the symbol '%s', but it does not exist" \
                "in this registry." % symbol)

        if hasattr(base_value, "in_base"):
            new_dimensions = base_value.units.dimensions
            base_value = base_value.in_base('cgs-ampere')
            base_value = base_value.value
        else:
            new_dimensions = self.lut[symbol][1]

        self.lut[symbol] = ((float(base_value), new_dimensions) +
                            self.lut[symbol][2:])

        if symbol in self.unit_objs:
            del self.unit_objs[symbol]

    def keys(self):
        """
        Print out the units contained in the lookup table.

        """
        return self.lut.keys()

    def to_json(self):
        """
        Returns a json-serialized version of the unit registry
        """
        sanitized_lut = {}
        for k, v in six.iteritems(self.lut):
            san_v = list(v)
            repr_dims = srepr(v[1])
            if SYMPY_VERSION < LooseVersion("1.0.0"):
                # see https://github.com/sympy/sympy/issues/6131
                repr_dims = re.sub("Symbol\('\([a-z_]*\)'\)",
                                   positive_symbol_replacer, repr_dims)
            san_v[1] = repr_dims
            sanitized_lut[k] = tuple(san_v)

        return json.dumps(sanitized_lut)

    @classmethod
    def from_json(cls, json_text):
        """
        Returns a UnitRegistry object from a json-serialized unit registry
        """
        data = json.loads(json_text)
        lut = {}
        for k, v in six.iteritems(data):
            unsan_v = list(v)
            unsan_v[1] = sympify(v[1])
            lut[k] = tuple(unsan_v)

        return cls(lut=lut, add_default_symbols=False)

    def list_same_dimensions(self, unit_object):
        """
        Return a list of base unit names that this registry knows about that
        are of equivalent dimensions to *unit_object*.
        """
        equiv = [k for k, v in self.lut.items()
                 if v[1] is unit_object.dimensions]
        equiv += [n for n, u in self.unit_objs.items()
                 if u.dimensions is unit_object.dimensions]
        equiv = list(sorted(set(equiv)))
        return equiv
