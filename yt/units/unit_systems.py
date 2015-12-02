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
from yt.units import dimensions
from yt.units.unit_object import Unit

class UnitSystem(object):
    def __init__(self, length_unit, mass_unit, time_unit,
                 temperature_unit, angle_unit):
        self.units_map = {dimensions.length: length_unit,
                          dimensions.mass: mass_unit,
                          dimensions.time: time_unit,
                          dimensions.temperature: temperature_unit,
                          dimensions.angle: angle_unit}

    def __getitem__(self, key):
        if isinstance(key, string_types):
            key = getattr(dimensions, key)
        if key not in self.units_map:
            dims = key.expand()
            units = Unit("")
            for factor in dims.as_ordered_factors():
                dim = list(factor.free_symbols)[0]
                u = self.units_map[dim]
                if factor.is_Pow:
                    u = u ** factor.as_base_exp()[1]
                units *= u
            self.units_map[key] = units
        return self.units_map[key]

    def __setitem__(self, key, value):
        if isinstance(key, string_types):
            key = getattr(dimensions, key)
        self.units_map[key] = Unit(value)

cgs_unit_system = UnitSystem("cm", "g", "s", "K", "radian")
cgs_unit_system[dimensions.energy] = "erg"
cgs_unit_system[dimensions.pressure] = "dyne/cm**2"
cgs_unit_system[dimensions.force] = "dyne"
cgs_unit_system[dimensions.magnetic_field_cgs] = "gauss"

mks_unit_system = UnitSystem("m", "kg", "s", "K", "radian")
mks_unit_system[dimensions.current_mks] = "A"
mks_unit_system[dimensions.energy] = "J"
mks_unit_system[dimensions.pressure] = "Pa"
mks_unit_system[dimensions.force] = "N"
mks_unit_system[dimensions.magnetic_field_mks] = "T"

imperial_unit_system = UnitSystem("ft", "lbm", "s", "R", "radian")
imperial_unit_system[dimensions.force] = "lbf"

galactic_unit_system = UnitSystem("kpc", "Msun", "Myr", "K", "radian")
galactic_unit_system[dimensions.energy] = "keV"
galactic_unit_system[dimensions.magnetic_field_cgs] = "uG"

geographic_unit_system = UnitSystem("mile", "lbm", "s", "K", "radian")
geographic_unit_system[dimensions.pressure] = "atm"

unit_system_registry = {"cgs": cgs_unit_system,
                        "mks": mks_unit_system,
                        "imperial": imperial_unit_system,
                        "galactic": galactic_unit_system,
                        "geographic": geographic_unit_system}
