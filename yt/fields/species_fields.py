import re

import numpy as np

from yt.frontends.sph.data_structures import ParticleDataset
from yt.utilities.chemical_formulas import ChemicalFormula
from yt.utilities.physical_ratios import _primordial_mass_fraction

from .field_plugin_registry import register_field_plugin

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
        return data[ftype, f"{species}_density"] / data[ftype, "density"]

    return _frac


def _mass_from_cell_volume_and_density(ftype, species):
    def _mass(field, data):
        return data[ftype, f"{species}_density"] * data["index", "cell_volume"]

    return _mass


def _mass_from_particle_mass_and_fraction(ftype, species):
    def _mass(field, data):
        return data[ftype, f"{species}_fraction"] * data[ftype, "particle_mass"]

    return _mass


def _create_number_density_func(ftype, species):
    formula = ChemicalFormula(species)

    def _number_density(field, data):
        weight = formula.weight  # This is in AMU
        weight *= data.ds.units.physical_constants.amu_cgs
        return data[ftype, f"{species}_density"] / weight

    return _number_density


def _create_density_func(ftype, species):
    def _density(field, data):
        return data[ftype, f"{species}_fraction"] * data[ftype, "density"]

    return _density


def add_species_field_by_density(registry, ftype, species):
    """
    This takes a field registry, a fluid type, and a species name and then
    adds the other fluids based on that.  This assumes that the field
    "SPECIES_density" already exists and refers to mass density.
    """
    unit_system = registry.ds.unit_system

    registry.add_field(
        (ftype, f"{species}_fraction"),
        sampling_type="local",
        function=_create_fraction_func(ftype, species),
        units="",
    )

    if isinstance(registry.ds, ParticleDataset):
        _create_mass_func = _mass_from_particle_mass_and_fraction
    else:
        _create_mass_func = _mass_from_cell_volume_and_density
    registry.add_field(
        (ftype, f"{species}_mass"),
        sampling_type="local",
        function=_create_mass_func(ftype, species),
        units=unit_system["mass"],
    )

    registry.add_field(
        (ftype, f"{species}_number_density"),
        sampling_type="local",
        function=_create_number_density_func(ftype, species),
        units=unit_system["number_density"],
    )

    return [
        (ftype, f"{species}_number_density"),
        (ftype, f"{species}_density"),
        (ftype, f"{species}_mass"),
    ]


def add_species_field_by_fraction(registry, ftype, species):
    """
    This takes a field registry, a fluid type, and a species name and then
    adds the other fluids based on that.  This assumes that the field
    "SPECIES_fraction" already exists and refers to mass fraction.
    """
    unit_system = registry.ds.unit_system

    registry.add_field(
        (ftype, f"{species}_density"),
        sampling_type="local",
        function=_create_density_func(ftype, species),
        units=unit_system["density"],
    )

    if isinstance(registry.ds, ParticleDataset):
        _create_mass_func = _mass_from_particle_mass_and_fraction
    else:
        _create_mass_func = _mass_from_cell_volume_and_density
    registry.add_field(
        (ftype, f"{species}_mass"),
        sampling_type="local",
        function=_create_mass_func(ftype, species),
        units=unit_system["mass"],
    )

    registry.add_field(
        (ftype, f"{species}_number_density"),
        sampling_type="local",
        function=_create_number_density_func(ftype, species),
        units=unit_system["number_density"],
    )

    return [
        (ftype, f"{species}_number_density"),
        (ftype, f"{species}_density"),
        (ftype, f"{species}_mass"),
    ]


def add_species_aliases(registry, ftype, alias_species, species):
    r"""
    This takes a field registry, a fluid type, and two species names.
    The first species name is one you wish to alias to an existing species
    name.  For instance you might alias all "H_p0" fields to "H\_" fields
    to indicate that "H\_" fields are really just neutral hydrogen fields.
    This function registers field aliases for the density, number_density,
    mass, and fraction fields between the two species given in the arguments.
    """
    registry.alias((ftype, f"{alias_species}_density"), (ftype, f"{species}_density"))
    registry.alias((ftype, f"{alias_species}_fraction"), (ftype, f"{species}_fraction"))
    registry.alias(
        (ftype, f"{alias_species}_number_density"),
        (ftype, f"{species}_number_density"),
    )
    registry.alias((ftype, f"{alias_species}_mass"), (ftype, f"{species}_mass"))


def add_deprecated_species_aliases(registry, ftype, alias_species, species):
    """
    Add the species aliases but with deprecation warnings.
    """

    for suffix in ["density", "fraction", "number_density", "mass"]:
        add_deprecated_species_alias(registry, ftype, alias_species, species, suffix)


def add_deprecated_species_alias(registry, ftype, alias_species, species, suffix):
    """
    Add a deprecated species alias field.
    """

    unit_system = registry.ds.unit_system
    if suffix == "fraction":
        my_units = ""
    else:
        my_units = unit_system[suffix]

    def _dep_field(field, data):
        return data[ftype, f"{species}_{suffix}"]

    registry.add_field(
        (ftype, f"{alias_species}_{suffix}"),
        sampling_type="local",
        function=_dep_field,
        units=my_units,
    )


def add_nuclei_density_fields(registry, ftype):
    unit_system = registry.ds.unit_system
    elements = _get_all_elements(registry.species_names)
    for element in elements:
        registry.add_field(
            (ftype, f"{element}_nuclei_density"),
            sampling_type="local",
            function=_nuclei_density,
            units=unit_system["number_density"],
        )

    # Here, we add default nuclei and number density fields for H and
    # He if they are not defined above, and if it was requested by
    # setting "default_species_fields"
    if registry.ds.default_species_fields is None:
        return

    dsf = registry.ds.default_species_fields
    # Right now, this only handles default fields for H and He
    for element in ["H", "He"]:
        # If these elements are already present in the dataset,
        # DO NOT set them
        if element in elements:
            continue
        # First add the default nuclei density fields
        registry.add_field(
            (ftype, f"{element}_nuclei_density"),
            sampling_type="local",
            function=_default_nuclei_density,
            units=unit_system["number_density"],
        )
        # Set up number density fields for hydrogen, either fully ionized or neutral.
        if element == "H":
            if dsf == "ionized":
                state = "p1"
            elif dsf == "neutral":
                state = "p0"
            else:
                raise NotImplementedError(
                    f"'default_species_fields' option '{dsf}' is not implemented!"
                )
            registry.alias(
                (ftype, f"H_{state}_number_density"), (ftype, "H_nuclei_density")
            )
        # Set up number density fields for helium, either fully ionized or neutral.
        if element == "He":
            if dsf == "ionized":
                state = "p2"
            elif dsf == "neutral":
                state = "p0"
            registry.alias(
                (ftype, f"He_{state}_number_density"), (ftype, "He_nuclei_density")
            )
    # If we're fully ionized, we need to setup the electron number density field
    if (ftype, "El_number_density") not in registry and dsf == "ionized":
        registry.add_field(
            (ftype, "El_number_density"),
            sampling_type="local",
            function=_default_nuclei_density,
            units=unit_system["number_density"],
        )


def _default_nuclei_density(field, data):
    ftype = field.name[0]
    element = field.name[1][: field.name[1].find("_")]
    amu_cgs = data.ds.units.physical_constants.amu_cgs
    if element == "El":
        # This is for determining the electron number density.
        # If we got here, this assumes full ionization!
        muinv = 1.0 * _primordial_mass_fraction["H"] / ChemicalFormula("H").weight
        muinv += 2.0 * _primordial_mass_fraction["He"] / ChemicalFormula("He").weight
    else:
        # This is for anything else besides electrons
        muinv = _primordial_mass_fraction[element] / ChemicalFormula(element).weight
    return data[ftype, "density"] * muinv / amu_cgs


def _nuclei_density(field, data):
    ftype = field.name[0]
    element = field.name[1][: field.name[1].find("_")]

    nuclei_mass_field = f"{element}_nuclei_mass_density"
    if (ftype, nuclei_mass_field) in data.ds.field_info:
        return (
            data[(ftype, nuclei_mass_field)]
            / ChemicalFormula(element).weight
            / data.ds.units.physical_constants.amu_cgs
        )
    metal_field = f"{element}_metallicity"
    if (ftype, metal_field) in data.ds.field_info:
        return (
            data[ftype, "density"]
            * data[(ftype, metal_field)]
            / ChemicalFormula(element).weight
            / data.ds.units.physical_constants.amu_cgs
        )

    field_data = np.zeros_like(
        data[ftype, f"{data.ds.field_info.species_names[0]}_number_density"]
    )
    for species in data.ds.field_info.species_names:
        nucleus = species
        if "_" in species:
            nucleus = species[: species.find("_")]
        # num is the number of nuclei contributed by this species.
        num = _get_element_multiple(nucleus, element)
        # Since this is a loop over all species existing in this dataset,
        # we will encounter species that contribute nothing, so we skip them.
        if num == 0:
            continue
        field_data += num * data[ftype, f"{species}_number_density"]
    return field_data


def _get_all_elements(species_list):
    elements = []
    for species in species_list:
        for item in re.findall("[A-Z][a-z]?|[0-9]+", species):
            if not item.isdigit() and item not in elements and item != "El":
                elements.append(item)
    return elements


def _get_element_multiple(compound, element):
    my_split = re.findall("[A-Z][a-z]?|[0-9]+", compound)
    if element not in my_split:
        return 0
    loc = my_split.index(element)
    if loc == len(my_split) - 1 or not my_split[loc + 1].isdigit():
        return 1
    return int(my_split[loc + 1])


@register_field_plugin
def setup_species_fields(registry, ftype="gas", slice_info=None):
    for species in registry.species_names:
        # These are all the species we should be looking for fractions or
        # densities of.
        if (ftype, f"{species}_density") in registry:
            func = add_species_field_by_density
        elif (ftype, f"{species}_fraction") in registry:
            func = add_species_field_by_fraction
        else:
            # Skip it
            continue
        func(registry, ftype, species)

        # Add aliases of X_p0_<field> to X_<field>.
        # These are deprecated and will be removed soon.
        if ChemicalFormula(species).charge == 0:
            alias_species = species.split("_")[0]
            if (ftype, f"{alias_species}_density") in registry:
                continue
            add_deprecated_species_aliases(registry, "gas", alias_species, species)

    add_nuclei_density_fields(registry, ftype)
