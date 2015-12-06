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
from yt.units.unit_object import Unit, unit_system_registry

class UnitSystem(object):
    def __init__(self, name, length_unit, mass_unit, time_unit,
                 temperature_unit, angle_unit, current_mks_unit=None,
                 registry=None):
        self.registry = registry
        self.units_map = {dimensions.length: Unit(length_unit, registry=self.registry),
                          dimensions.mass: Unit(mass_unit, registry=self.registry),
                          dimensions.time: Unit(time_unit, registry=self.registry),
                          dimensions.temperature: Unit(temperature_unit, registry=self.registry),
                          dimensions.angle: Unit(angle_unit, registry=self.registry)}
        self._dims = ["length","mass","time","temperature","angle"]
        if current_mks_unit is not None:
            self.units_map[dimensions.current_mks] = Unit(current_mks_unit, registry=self.registry)
            self._dims.append("current_mks")
        self.registry = registry
        self.base_units = self.units_map.copy()
        unit_system_registry[name] = self
        self.name = name

    def __getitem__(self, key):
        if isinstance(key, string_types):
            if key not in self._dims:
                self._dims.append(key)
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
            if key not in self._dims:
                self._dims.append(key)
            key = getattr(dimensions, key)
        self.units_map[key] = Unit(value, registry=self.registry)

    def __str__(self):
        return self.name

    def __repr__(self):
        repr = "%s Unit System\n" % self.name
        repr += " Base Units:\n"
        for dim in self.base_units:
            repr += "  %s: %s\n" % (str(dim).strip("()"), self.base_units[dim])
        repr += " Other Units:\n"
        for key in self._dims:
            dim = getattr(dimensions, key)
            if dim not in self.base_units:
                repr += "  %s: %s\n" % (key, self.units_map[dim])
        return repr

def create_code_unit_system(ds):
    code_unit_system = UnitSystem(str(ds), "code_length", "code_mass", "code_time",
                                  "code_temperature", "radian", registry=ds.unit_registry)
    code_unit_system["velocity"] = "code_velocity"

cgs_unit_system = UnitSystem("cgs", "cm", "g", "s", "K", "radian")
cgs_unit_system["energy"] = "erg"
cgs_unit_system["specific_energy"] = "erg/g"
cgs_unit_system["pressure"] = "dyne/cm**2"
cgs_unit_system["force"] = "dyne"
cgs_unit_system["magnetic_field_cgs"] = "gauss"

mks_unit_system = UnitSystem("mks", "m", "kg", "s", "K", "radian",
                             current_mks_unit="A")
mks_unit_system["energy"] = "J"
mks_unit_system["specific_energy"] = "J/kg"
mks_unit_system["pressure"] = "Pa"
mks_unit_system["force"] = "N"
mks_unit_system["magnetic_field_mks"] = "T"

imperial_unit_system = UnitSystem("imperial", "ft", "lbm", "s", "R", "radian")
imperial_unit_system["force"] = "lbf"
imperial_unit_system["energy"] = "ft*lbf"
imperial_unit_system["pressure"] = "lbf/ft**2"

galactic_unit_system = UnitSystem("galactic", "kpc", "Msun", "Myr", "K", "radian")
galactic_unit_system["energy"] = "keV"
galactic_unit_system["magnetic_field_cgs"] = "uG"