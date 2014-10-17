"""
Base dimensions


"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import yt.utilities.physical_constants as pc
from yt.units.dimensions import temperature, mass, energy, length, rate
from yt.extern.six import add_metaclass

equivalence_registry = {}

class RegisteredEquivalence(type):
    def __init__(cls, name, b, d):
        type.__init__(cls, name, b, d)
        if hasattr(cls, "_type_name") and not cls._skip_add:
            equivalence_registry[cls._type_name] = cls

@add_metaclass(RegisteredEquivalence)
class Equivalence(object):
    _skip_add = False

class ThermalEquivalence(Equivalence):
    _type_name = "thermal"
    dims = (temperature,energy,)

    def convert(self, x, new_dims):
        if new_dims == energy:
            return pc.kboltz*x
        elif new_dims == temperature:
            return x/pc.kboltz

    def __str__(self):
        return "thermal: temperature <-> energy"

class MassEnergyEquivalence(Equivalence):
    _type_name = "mass_energy"
    dims = (mass,energy,)

    def convert(self, x, new_dims):
        if new_dims == energy:
            return x*pc.clight*pc.clight
        elif new_dims == mass:
            return x/(pc.clight*pc.clight)

    def __str__(self):
        return "mass_energy: mass <-> energy"

class SpectralEquivalence(Equivalence):
    _type_name = "spectral"
    dims = (length,rate,energy,)

    def convert(self, x, new_dims):
        if new_dims == energy:
            if x.units.dimensions == length:
                nu = pc.clight/x
            elif x.units.dimensions == rate:
                nu = x
            return pc.hcgs*nu
        elif new_dims == length:
            if x.units.dimensions == rate:
                return pc.clight/x
            elif x.units.dimensions == energy:
                return pc.hcgs*pc.clight/x
        elif new_dims == rate:
            if x.units.dimensions == length:
                return pc.clight/x
            elif x.units.dimensions == energy:
                return x/pc.hcgs

    def __str__(self):
        return "spectral: length <-> rate <-> energy"

