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
import re

from yt.utilities.physical_constants import \
    amu_cgs
from yt.utilities.physical_ratios import \
    primordial_H_mass_fraction

from yt.utilities.chemical_formulas import \
    ChemicalFormula
from .field_plugin_registry import \
    register_field_plugin

_primordial_mass_fraction = \
  {"H": primordial_H_mass_fraction,
   "He" : (1 - primordial_H_mass_fraction)}

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
             / weight
    return _number_density

def _create_density_func(ftype, species):
    def _density(field, data):
        return data[ftype, "%s_fraction" % species] \
            * data[ftype,'density']
    return _density

def add_species_field_by_density(registry, ftype, species, 
                                 particle_type = False):
    """
    This takes a field registry, a fluid type, and a species name and then
    adds the other fluids based on that.  This assumes that the field
    "SPECIES_density" already exists and refers to mass density.
    """
    unit_system = registry.ds.unit_system

    registry.add_field((ftype, "%s_fraction" % species), 
                       function = _create_fraction_func(ftype, species),
                       particle_type = particle_type,
                       units = "")

    registry.add_field((ftype, "%s_mass" % species),
                       function = _create_mass_func(ftype, species),
                       particle_type = particle_type,
                       units = unit_system["mass"])

    registry.add_field((ftype, "%s_number_density" % species),
                       function = _create_number_density_func(ftype, species),
                       particle_type = particle_type,
                       units = unit_system["number_density"])

    return [(ftype, "%s_number_density" % species),
            (ftype, "%s_density" % species),
            (ftype, "%s_mass" % species)]

def add_species_field_by_fraction(registry, ftype, species, 
                                  particle_type = False):
    """
    This takes a field registry, a fluid type, and a species name and then
    adds the other fluids based on that.  This assumes that the field
    "SPECIES_fraction" already exists and refers to mass fraction.
    """
    unit_system = registry.ds.unit_system

    registry.add_field((ftype, "%s_density" % species), 
                       function = _create_density_func(ftype, species),
                       particle_type = particle_type,
                       units = unit_system["density"])

    registry.add_field((ftype, "%s_mass" % species),
                       function = _create_mass_func(ftype, species),
                       particle_type = particle_type,
                       units = unit_system["mass"])

    registry.add_field((ftype, "%s_number_density" % species),
                       function = _create_number_density_func(ftype, species),
                       particle_type = particle_type,
                       units = unit_system["number_density"])

    return [(ftype, "%s_number_density" % species),
            (ftype, "%s_density" % species),
            (ftype, "%s_mass" % species)]

def add_species_aliases(registry, ftype, alias_species, species):
    """
    This takes a field registry, a fluid type, and two species names.  
    The first species name is one you wish to alias to an existing species
    name.  For instance you might alias all "H_p0" fields to "H_" fields
    to indicate that "H_" fields are really just neutral hydrogen fields.
    This function registers field aliases for the density, number_density,
    mass, and fraction fields between the two species given in the arguments.
    """
    registry.alias((ftype, "%s_density" % alias_species), 
                   (ftype, "%s_density" % species))
    registry.alias((ftype, "%s_fraction" % alias_species), 
                   (ftype, "%s_fraction" % species))
    registry.alias((ftype, "%s_number_density" % alias_species), 
                   (ftype, "%s_number_density" % species))
    registry.alias((ftype, "%s_mass" % alias_species), 
                   (ftype, "%s_mass" % species))

def add_nuclei_density_fields(registry, ftype,
                              particle_type = False):
    unit_system = registry.ds.unit_system
    elements = _get_all_elements(registry.species_names)
    for element in elements:
        registry.add_field((ftype, "%s_nuclei_density" % element),
                           function = _nuclei_density,
                           particle_type = particle_type,
                           units = unit_system["number_density"])

    for element in ["H", "He"]:
        if element in elements:
            continue
        registry.add_field((ftype, "%s_nuclei_density" % element),
                           function = _default_nuclei_density,
                           particle_type = particle_type,
                           units = unit_system["number_density"])

def _default_nuclei_density(field, data):
    ftype = field.name[0]
    element = field.name[1][:field.name[1].find("_")]
    return data[ftype, "density"] * _primordial_mass_fraction[element] / \
      ChemicalFormula(element).weight / amu_cgs
        
def _nuclei_density(field, data):
    ftype = field.name[0]
    element = field.name[1][:field.name[1].find("_")]

    nuclei_mass_field = "%s_nuclei_mass_density" % element
    if (ftype, nuclei_mass_field) in data.ds.field_info:
        return data[(ftype, nuclei_mass_field)] / \
          ChemicalFormula(element).weight / amu_cgs
    metal_field = "%s_metallicity" % element
    if (ftype, metal_field) in data.ds.field_info:
        return data[ftype, "density"] * data[(ftype, metal_field)] / \
          ChemicalFormula(element).weight / amu_cgs

    field_data = np.zeros_like(data[ftype, "%s_number_density" %
                                    data.ds.field_info.species_names[0]])
    for species in data.ds.field_info.species_names:
        nucleus = species
        if "_" in species:
            nucleus = species[:species.find("_")]
        # num is the number of nuclei contributed by this species.
        num = _get_element_multiple(nucleus, element)
        # Since this is a loop over all species existing in this dataset,
        # we will encounter species that contribute nothing, so we skip them.
        if num == 0:
            continue
        field_data += num * data[ftype, "%s_number_density" % species]
    return field_data

def _get_all_elements(species_list):
    elements = []
    for species in species_list:
        for item in re.findall('[A-Z][a-z]?|[0-9]+', species):
            if not item.isdigit() and item not in elements \
              and item != "El":
                elements.append(item)
    return elements
    
def _get_element_multiple(compound, element):
    my_split = re.findall('[A-Z][a-z]?|[0-9]+', compound)
    if element not in my_split:
        return 0
    loc = my_split.index(element)
    if loc == len(my_split) - 1 or not my_split[loc + 1].isdigit():
        return 1
    return int(my_split[loc + 1])

@register_field_plugin
def setup_species_fields(registry, ftype = "gas", slice_info = None):
    # We have to check what type of field this is -- if it's particles, then we
    # set particle_type to True.
    particle_type = ftype not in registry.ds.fluid_types
    for species in registry.species_names:
        # These are all the species we should be looking for fractions or
        # densities of.
        if (ftype, "%s_density" % species) in registry:
            func = add_species_field_by_density
        elif (ftype, "%s_fraction" % species) in registry:
            func = add_species_field_by_fraction
        else:
            # Skip it
            continue
        func(registry, ftype, species, particle_type)
        # Adds aliases for all neutral species from their raw "MM_"
        # species to "MM_p0_" species to be explicit.
        # See YTEP-0003 for more details.
        if (ChemicalFormula(species).charge == 0):
            alias_species = "%s_p0" % species.split('_')[0]
            add_species_aliases(registry, "gas", alias_species, species)
    add_nuclei_density_fields(registry, ftype, particle_type=particle_type)
