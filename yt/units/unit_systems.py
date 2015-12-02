"""
Unit system class.

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.extern.six import string_types
from yt.units.unit_lookup_table import cgs_base_units, mks_base_units
from yt.units import dimensions
from yt.units.unit_object import Unit

class UnitSystem(object):
    def __init__(self, base_units):
        self.base_units = {}
        for key, value in base_units.items():
            self.base_units[key] = Unit(value)

    def __getitem__(self, key):
        if isinstance(key, string_types):
            key = getattr(dimensions, key)
        if key not in self.base_units:
            dims = key.expand()
            units = Unit("")
            for factor in dims.as_ordered_factors():
                dim = list(factor.free_symbols)[0]
                u = self.base_units[dim]
                if factor.is_Pow:
                    u = u ** factor.as_base_exp()[1]
                units *= u
            self.base_units[key] = units
        return self.base_units[key]

    def __setitem__(self, key, value):
        if isinstance(key, string_types):
            key = getattr(dimensions, key)
        self.base_units[key] = Unit(value)

cgs_unit_system = UnitSystem(cgs_base_units)
cgs_unit_system[dimensions.energy] = "erg"
cgs_unit_system[dimensions.pressure] = "dyne/cm**2"
cgs_unit_system[dimensions.force] = "dyne"
cgs_unit_system[dimensions.magnetic_field_cgs] = "gauss"

mks_unit_system = UnitSystem(mks_base_units)
mks_unit_system[dimensions.energy] = "J"
mks_unit_system[dimensions.pressure] = "Pa"
mks_unit_system[dimensions.force] = "N"
mks_unit_system[dimensions.magnetic_field_mks] = "T"