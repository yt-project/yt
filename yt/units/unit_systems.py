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
from yt.units.unit_object import Unit, unit_system_registry, _get_system_unit_string
from yt.utilities import physical_constants as pc


class UnitSystemConstants(object):
    """
    A class to facilitate conversions of physical constants into a given unit
    system specified by *name*.
    """
    def __init__(self, name):
        self.name = name

    def __repr__(self):
        return "Physical constants in %s units." % self.name

    def __str__(self):
        return self.name

    def __getattr__(self, item):
        return getattr(pc, item).in_base(self.name)


class UnitSystem(object):
    """
    Create a UnitSystem for facilitating conversions to a default set of units.

    Parameters
    ----------
    name : string
        The name of the unit system. Will be used as the key in the 
        *unit_system_registry* dict to reference the unit system by.
    length_unit : string
        The base length unit of this unit system.
    mass_unit : string
        The base mass unit of this unit system.
    time_unit : string
        The base time unit of this unit system.
    temperature_unit : string, optional
        The base temperature unit of this unit system. Defaults to "K".
    angle_unit : string, optional
        The base angle unit of this unit system. Defaults to "rad".
    current_mks_unit : string, optional
        The base current unit of this unit system. Only used in MKS 
        or MKS-based unit systems.
    registry : :class:`yt.units.unit_registry.UnitRegistry` object
        The unit registry associated with this unit system. Only 
        useful for defining unit systems based on code units.
    """
    def __init__(self, name, length_unit, mass_unit, time_unit,
                 temperature_unit="K", angle_unit="rad", current_mks_unit=None,
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
        self.constants = UnitSystemConstants(self.name)

    def __getitem__(self, key):
        if isinstance(key, string_types):
            key = getattr(dimensions, key)
        um = self.units_map
        if key not in um or um[key].dimensions is not key:
            units = _get_system_unit_string(key, self.units_map)
            self.units_map[key] = Unit(units, registry=self.registry)
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


def _make_unit_system_copy(unit_registry, unit_system):
    """
    Make a copy of a unit system with a different unit registry.

    Parameters
    ----------
    unit_registry : UnitRegistry instance
        The unit registry we want to use with this unit system,
        most likely from a dataset.
    unit_system : string
        The name of the unit system we want to make a copy of 
        with the new registry.
    """
    unit_system = unit_system_registry[unit_system]
    base_units = ["length", "mass", "time", "temperature", "angle"]
    if "current_mks" in unit_system._dims:
        current_mks_unit = str(unit_system["current_mks"])
        base_units.append("current_mks")
    else:
        current_mks_unit = None
    sys_name = "{}_{}".format(unit_system.name, unit_registry.unit_system_id)
    ds_unit_system = UnitSystem(sys_name, str(unit_system["length"]), 
                                str(unit_system["mass"]), 
                                str(unit_system["time"]), 
                                temperature_unit=str(unit_system["temperature"]),
                                angle_unit=str(unit_system["angle"]),
                                current_mks_unit=current_mks_unit,
                                registry=unit_registry)
    for dim in unit_system._dims:
        if dim not in base_units:
            ds_unit_system[dim] = str(unit_system[dim])
    return ds_unit_system


def create_code_unit_system(unit_registry, current_mks_unit=None):
    code_unit_system = UnitSystem(unit_registry.unit_system_id, "code_length", 
                                  "code_mass", "code_time", "code_temperature",
                                  current_mks_unit=current_mks_unit, 
                                  registry=unit_registry)
    code_unit_system["velocity"] = "code_velocity"
    if current_mks_unit:
        code_unit_system["magnetic_field_mks"] = "code_magnetic"
    else:
        code_unit_system["magnetic_field_cgs"] = "code_magnetic"
    code_unit_system["pressure"] = "code_pressure"


cgs_unit_system = UnitSystem("cgs", "cm", "g", "s")
cgs_unit_system["energy"] = "erg"
cgs_unit_system["specific_energy"] = "erg/g"
cgs_unit_system["pressure"] = "dyne/cm**2"
cgs_unit_system["force"] = "dyne"
cgs_unit_system["magnetic_field_cgs"] = "gauss"
cgs_unit_system["charge_cgs"] = "esu"
cgs_unit_system["current_cgs"] = "statA"

mks_unit_system = UnitSystem("mks", "m", "kg", "s", current_mks_unit="A")
mks_unit_system["energy"] = "J"
mks_unit_system["specific_energy"] = "J/kg"
mks_unit_system["pressure"] = "Pa"
mks_unit_system["force"] = "N"
mks_unit_system["magnetic_field_mks"] = "T"
mks_unit_system["charge_mks"] = "C"

imperial_unit_system = UnitSystem("imperial", "ft", "lbm", "s", 
                                  temperature_unit="R")
imperial_unit_system["force"] = "lbf"
imperial_unit_system["energy"] = "ft*lbf"
imperial_unit_system["pressure"] = "lbf/ft**2"

galactic_unit_system = UnitSystem("galactic", "kpc", "Msun", "Myr")
galactic_unit_system["energy"] = "keV"
galactic_unit_system["magnetic_field_cgs"] = "uG"

solar_unit_system = UnitSystem("solar", "AU", "Mearth", "yr")

geometrized_unit_system = UnitSystem("geometrized", "l_geom", 
                                     "m_geom", "t_geom")

planck_unit_system = UnitSystem("planck", "l_pl", "m_pl", "t_pl", 
                                temperature_unit="T_pl")
planck_unit_system["energy"] = "E_pl"
planck_unit_system["charge_cgs"] = "q_pl"


cgs_ampere_unit_system = UnitSystem('cgs-ampere', 'cm', 'g', 's',
                                    current_mks_unit='A')
cgs_ampere_unit_system["energy"] = "erg"
cgs_ampere_unit_system["specific_energy"] = "erg/g"
cgs_ampere_unit_system["pressure"] = "dyne/cm**2"
cgs_ampere_unit_system["force"] = "dyne"
cgs_ampere_unit_system["magnetic_field_cgs"] = "gauss"
cgs_ampere_unit_system["charge_cgs"] = "esu"
cgs_ampere_unit_system["current_cgs"] = "statA"
