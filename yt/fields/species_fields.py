import numpy as np
import re

from yt.fields.field_detector import \
    FieldDetector
from yt.frontends.sph.data_structures import \
    ParticleDataset
from yt.funcs import \
    issue_deprecation_warning
from yt.utilities.physical_ratios import \
    _primordial_mass_fraction
from yt.utilities.chemical_formulas import \
    ChemicalFormula
from .field_plugin_registry import \
    register_field_plugin


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

def _mass_from_cell_volume_and_density(ftype, species):
    def _mass(field, data):
        return data[ftype, "%s_density" % species] \
             * data["index", "cell_volume"]
    return _mass

def _mass_from_particle_mass_and_fraction(ftype, species):
    def _mass(field, data):
        return data[ftype, "%s_fraction" % species] \
            * data[ftype, 'particle_mass']
    return _mass

def _create_number_density_func(ftype, species):
    formula = ChemicalFormula(species)
    def _number_density(field, data):
        weight = formula.weight # This is in AMU
        weight *= data.ds.units.physical_constants.amu_cgs
        return data[ftype, "%s_density" % species] \
             / weight
    return _number_density

def _create_density_func(ftype, species):
    def _density(field, data):
        return data[ftype, "%s_fraction" % species] \
            * data[ftype,'density']
    return _density

def add_species_field_by_density(registry, ftype, species):
    """
    This takes a field registry, a fluid type, and a species name and then
    adds the other fluids based on that.  This assumes that the field
    "SPECIES_density" already exists and refers to mass density.
    """
    unit_system = registry.ds.unit_system

    registry.add_field((ftype, "%s_fraction" % species),
                       sampling_type="local",
                       function = _create_fraction_func(ftype, species),
                       units = "")

    if isinstance(registry.ds, ParticleDataset):
        _create_mass_func = _mass_from_particle_mass_and_fraction
    else:
        _create_mass_func = _mass_from_cell_volume_and_density
    registry.add_field((ftype, "%s_mass" % species),
                       sampling_type="local",
                       function = _create_mass_func(ftype, species),
                       units = unit_system["mass"])

    registry.add_field((ftype, "%s_number_density" % species),
                       sampling_type="local",
                       function = _create_number_density_func(ftype, species),
                       units = unit_system["number_density"])

    return [(ftype, "%s_number_density" % species),
            (ftype, "%s_density" % species),
            (ftype, "%s_mass" % species)]


def add_species_field_by_fraction(registry, ftype, species):
    """
    This takes a field registry, a fluid type, and a species name and then
    adds the other fluids based on that.  This assumes that the field
    "SPECIES_fraction" already exists and refers to mass fraction.
    """
    unit_system = registry.ds.unit_system

    registry.add_field((ftype, "%s_density" % species),
                       sampling_type="local",
                       function = _create_density_func(ftype, species),
                       units = unit_system["density"])

    if isinstance(registry.ds, ParticleDataset):
        _create_mass_func = _mass_from_particle_mass_and_fraction
    else:
        _create_mass_func = _mass_from_cell_volume_and_density
    registry.add_field((ftype, "%s_mass" % species),
                       sampling_type="local",
                       function = _create_mass_func(ftype, species),
                       units = unit_system["mass"])

    registry.add_field((ftype, "%s_number_density" % species),
                       sampling_type="local",
                       function = _create_number_density_func(ftype, species),
                       units = unit_system["number_density"])

    return [(ftype, "%s_number_density" % species),
            (ftype, "%s_density" % species),
            (ftype, "%s_mass" % species)]

def add_species_aliases(registry, ftype, alias_species, species):
    """
    This takes a field registry, a fluid type, and two species names.  
    The first species name is one you wish to alias to an existing species
    name.  For instance you might alias all "H_p0" fields to "H\_" fields
    to indicate that "H\_" fields are really just neutral hydrogen fields.
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

def add_deprecated_species_aliases(registry, ftype, alias_species, species):
    """
    Add the species aliases but with deprecation warnings.
    """

    for suffix in ["density", "fraction", "number_density", "mass"]:
        add_deprecated_species_alias(
            registry, ftype, alias_species, species, suffix)

def add_deprecated_species_alias(registry, ftype, alias_species, species,
                                 suffix):
    """
    Add a deprecated species alias field.
    """

    unit_system = registry.ds.unit_system
    if suffix == "fraction":
        my_units = ""
    else:
        my_units = unit_system[suffix]

    def _dep_field(field, data):
        if not isinstance(data, FieldDetector):
            issue_deprecation_warning(
                ("The \"%s_%s\" field is deprecated. " +
                 "Please use \"%s_%s\" instead.") %
                (alias_species, suffix, species, suffix))
        return data[ftype, "%s_%s" % (species, suffix)]

    registry.add_field((ftype, "%s_%s" % (alias_species, suffix)),
                       sampling_type="local",
                       function=_dep_field,
                       units=my_units)

def add_nuclei_density_fields(registry, ftype):
    unit_system = registry.ds.unit_system
    elements = _get_all_elements(registry.species_names)
    for element in elements:
        registry.add_field((ftype, "%s_nuclei_density" % element),
                           sampling_type="local",
                           function=_nuclei_density,
                           units=unit_system["number_density"])

    # Here, we add default nuclei and number density fields for H and
    # He if they are not defined above. This assumes full ionization!
    for element in ["H", "He"]:
        if element in elements:
            continue
        registry.add_field((ftype, "%s_nuclei_density" % element), 
                           sampling_type="local",
                           function=_default_nuclei_density,
                           units=unit_system["number_density"])
        if element == "H":
            registry.alias((ftype, "H_p1_number_density"),
                           (ftype, "H_nuclei_density"))

        if element == "He":
            registry.alias((ftype, "He_p2_number_density"),
                           (ftype, "He_nuclei_density"))

    if (ftype, "El_number_density") not in registry:
        registry.add_field((ftype, "El_number_density"),
                           sampling_type="local",
                           function=_default_nuclei_density,
                           units=unit_system["number_density"])


def _default_nuclei_density(field, data):
    ftype = field.name[0]
    element = field.name[1][:field.name[1].find("_")]
    amu_cgs = data.ds.units.physical_constants.amu_cgs
    if element == "El":
        # This assumes full ionization!
        muinv = 1.0*_primordial_mass_fraction["H"] / \
          ChemicalFormula("H").weight
        muinv += 2.0*_primordial_mass_fraction["He"] / \
          ChemicalFormula("He").weight
    else:
        muinv = _primordial_mass_fraction[element] / \
          ChemicalFormula(element).weight
    return data[ftype, "density"] * muinv / amu_cgs


def _nuclei_density(field, data):
    ftype = field.name[0]
    element = field.name[1][:field.name[1].find("_")]

    nuclei_mass_field = "%s_nuclei_mass_density" % element
    if (ftype, nuclei_mass_field) in data.ds.field_info:
        return data[(ftype, nuclei_mass_field)] / \
          ChemicalFormula(element).weight / data.ds.units.physical_constants.amu_cgs
    metal_field = "%s_metallicity" % element
    if (ftype, metal_field) in data.ds.field_info:
        return data[ftype, "density"] * data[(ftype, metal_field)] / \
          ChemicalFormula(element).weight / data.ds.units.physical_constants.amu_cgs

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
        func(registry, ftype, species)

        # Add aliases of X_p0_<field> to X_<field>.
        # These are deprecated and will be removed soon.
        if ChemicalFormula(species).charge == 0:
            alias_species = species.split("_")[0]
            if (ftype, "{}_density".format(alias_species)) in registry:
                continue
            add_deprecated_species_aliases(
                registry, "gas", alias_species, species)

    add_nuclei_density_fields(registry, ftype)
