"""
Equivalencies between different kinds of units

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.units.dimensions import temperature, mass, energy, length, rate, \
    velocity, dimensionless, density, number_density, flux, current_cgs, \
    current_mks, charge_cgs, charge_mks, magnetic_field_cgs, magnetic_field_mks, \
    electric_potential_cgs, electric_potential_mks, electric_field_cgs, \
    electric_field_mks, resistance_cgs, resistance_mks

from yt.utilities.physical_ratios import speed_of_light_cm_per_s
from yt.extern.six import add_metaclass
import numpy as np

equivalence_registry = {}

class RegisteredEquivalence(type):
    def __init__(cls, name, b, d):
        type.__init__(cls, name, b, d)
        if hasattr(cls, "_type_name") and not cls._skip_add:
            equivalence_registry[cls._type_name] = cls

@add_metaclass(RegisteredEquivalence)
class Equivalence(object):
    _skip_add = False
    _one_way = False
    def __init__(self):
        import yt.utilities.physical_constants as pc
        self.pc = pc

class NumberDensityEquivalence(Equivalence):
    _type_name = "number_density"
    dims = (density,number_density,)

    def convert(self, x, new_dims, mu=0.6):
        if new_dims == number_density:
            return x/(mu*self.pc.mh)
        elif new_dims == density:
            return x*mu*self.pc.mh

    def __str__(self):
        return "number density: density <-> number density"

class ThermalEquivalence(Equivalence):
    _type_name = "thermal"
    dims = (temperature,energy,)

    def convert(self, x, new_dims):
        if new_dims == energy:
            return self.pc.kboltz*x
        elif new_dims == temperature:
            return x/self.pc.kboltz

    def __str__(self):
        return "thermal: temperature <-> energy"

class MassEnergyEquivalence(Equivalence):
    _type_name = "mass_energy"
    dims = (mass,energy,)

    def convert(self, x, new_dims):
        if new_dims == energy:
            return x*self.pc.clight*self.pc.clight
        elif new_dims == mass:
            return x/(self.pc.clight*self.pc.clight)

    def __str__(self):
        return "mass_energy: mass <-> energy"

class SpectralEquivalence(Equivalence):
    _type_name = "spectral"
    dims = (length,rate,energy,)

    def convert(self, x, new_dims):
        if new_dims == energy:
            if x.units.dimensions == length:
                nu = self.pc.clight/x
            elif x.units.dimensions == rate:
                nu = x
            return self.pc.hcgs*nu
        elif new_dims == length:
            if x.units.dimensions == rate:
                return self.pc.clight/x
            elif x.units.dimensions == energy:
                return self.pc.hcgs*self.pc.clight/x
        elif new_dims == rate:
            if x.units.dimensions == length:
                return self.pc.clight/x
            elif x.units.dimensions == energy:
                return x/self.pc.hcgs

    def __str__(self):
        return "spectral: length <-> rate <-> energy"

class SoundSpeedEquivalence(Equivalence):
    _type_name = "sound_speed"
    dims = (velocity,temperature,energy,)

    def convert(self, x, new_dims, mu=0.6, gamma=5./3.):
        if new_dims == velocity:
            if x.units.dimensions == temperature:
                kT = self.pc.kboltz*x
            elif x.units.dimensions == energy:
                kT = x
            return np.sqrt(gamma*kT/(mu*self.pc.mh))
        else:
            kT = x*x*mu*self.pc.mh/gamma
            if new_dims == temperature:
                return kT/self.pc.kboltz
            else:
                return kT

    def __str__(self):
        return "sound_speed (ideal gas): velocity <-> temperature <-> energy"

class LorentzEquivalence(Equivalence):
    _type_name = "lorentz"
    dims = (dimensionless,velocity,)

    def convert(self, x, new_dims):
        if new_dims == dimensionless:
            beta = x.in_cgs()/self.pc.clight
            return 1./np.sqrt(1.-beta**2)
        elif new_dims == velocity:
            return self.pc.clight*np.sqrt(1.-1./(x*x))

    def __str__(self):
        return "lorentz: velocity <-> dimensionless"

class SchwarzschildEquivalence(Equivalence):
    _type_name = "schwarzschild"
    dims = (mass,length,)

    def convert(self, x, new_dims):
        if new_dims == length:
            return 2.*self.pc.G*x/(self.pc.clight*self.pc.clight)
        elif new_dims == mass:
            return 0.5*x*self.pc.clight*self.pc.clight/self.pc.G

    def __str__(self):
        return "schwarzschild: mass <-> length"

class ComptonEquivalence(Equivalence):
    _type_name = "compton"
    dims = (mass,length,)

    def convert(self, x, new_dims):
        return self.pc.hcgs/(x*self.pc.clight)

    def __str__(self):
        return "compton: mass <-> length"

class EffectiveTemperature(Equivalence):
    _type_name = "effective_temperature"
    dims = (flux,temperature,)

    def convert(self, x, new_dims):
        if new_dims == flux:
            return self.pc.stefan_boltzmann_constant_cgs*x**4
        elif new_dims == temperature:
            return (x/self.pc.stefan_boltzmann_constant_cgs)**0.25

    def __str__(self):
        return "effective_temperature: flux <-> temperature"

em_conversions = {
    charge_mks:("esu", 0.1*speed_of_light_cm_per_s),
    magnetic_field_mks:("gauss", 1.0e4),
    current_mks:("statA", 0.1*speed_of_light_cm_per_s),
    electric_potential_mks:("statV", 1.0e-8*speed_of_light_cm_per_s),
    resistance_mks:("statohm", 1.0e9/(speed_of_light_cm_per_s**2)),
    charge_cgs:("C", 10.0/speed_of_light_cm_per_s),
    magnetic_field_cgs:("T", 1.0e-4),
    current_cgs:("A", 10.0/speed_of_light_cm_per_s),
    electric_potential_cgs:("V", 1.0e8/speed_of_light_cm_per_s),
    resistance_cgs:("ohm", speed_of_light_cm_per_s**2*1.0e-9),
}

class ElectromagneticSI(Equivalence):
    _type_name = "SI"
    _one_way = True
    dims = (current_cgs, charge_cgs, magnetic_field_cgs,
            electric_field_cgs, electric_potential_cgs,
            resistance_cgs)
    
    def convert(self, x, new_dims):
        old_dims = x.units.dimensions
        new_units, convert_factor = em_conversions[old_dims]
        return x.in_cgs().v*convert_factor, new_units
    
    def __str__(self):
        return "SI: EM CGS unit -> EM SI unit"

class ElectromagneticCGS(Equivalence):
    _type_name = "CGS"
    _one_way = True
    dims = (current_mks, charge_mks, magnetic_field_mks,
            electric_field_mks, electric_potential_mks,
            resistance_mks)

    def convert(self, x, new_dims):
        old_dims = x.units.dimensions
        new_units, convert_factor = em_conversions[old_dims]
        return x.in_mks().v*convert_factor, new_units

    def __str__(self):
        return "CGS: EM SI unit -> EM CGS unit"

