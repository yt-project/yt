"""
Fields based on species of molecules or atoms.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from yt.utilities.physical_constants import \
    mh, \
    mass_sun_cgs, \
    amu_cgs
from yt.funcs import *
from yt.utilities.chemical_formulas import \
    ChemicalFormula

# See YTEP-0003 for details, but we want to ensure these fields are all
# populated:
#
#   * _mass
#   * _density
#   * _fraction
#   * _number_density
#

def _create_fraction_func(ftype, species):
    def _frac(field, data):
        return data[ftype, "%s_density" % species] \
             / data[ftype, "density"]
    return _frac

def _create_mass_func(ftype, species):
    def _mass(field, data):
        return data[ftype, "%s_density" % species] \
             * data["index", "cell_volume"]
    return _mass

def _create_number_density_func(ftype, species):
    formula = ChemicalFormula(species)
    weight = formula.weight # This is in AMU
    weight *= amu_cgs
    def _number_density(field, data):
        return data[ftype, "%s_density" % species] \
             / amu_cgs
    return _number_density

def _create_density_function(ftype, species):
    def _density(field, data):
        return data[ftype, "%s_fraction" % species]
    return _density

def add_species_field_by_density(registry, ftype, species):
    """
    This takes a field registry, a fluid type, and a species name and then
    adds the other fluids based on that.  This assumes that the field
    "SPECIES_density" already exists and refers to mass density.
    """
    registry.add_field((ftype, "%s_fraction" % species), 
                        function = _create_fraction_func(ftype, species),
                        units = "")
    registry.add_field((ftype, "%s_mass" % species),
                        function = _create_mass_func(ftype, species),
                        units = "g")
    registry.add_field((ftype, "%s_number_density" % species),
                        function = _create_number_density_func(ftype, species),
                        units = "cm**-3")

def add_species_field_by_fraction(registry, ftype, species):
    """
    This takes a field registry, a fluid type, and a species name and then
    adds the other fluids based on that.  This assumes that the field
    "SPECIES_fraction" already exists and refers to mass fraction.
    """
    registry.add_field((ftype, "%s_density" % species), 
                        function = _create_density_func(ftype, species),
                        units = "g/cm**3")
    registry.add_field((ftype, "%s_mass" % species),
                        function = _create_mass_func(ftype, species),
                        units = "g")
    registry.add_field((ftype, "%s_number_density" % species),
                        function = _create_number_density_func(ftype, species),
                        units = "cm**-3")
