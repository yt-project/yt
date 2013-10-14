"""
Fields specific to Enzo



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from yt.data_objects.field_info_container import \
    FieldInfoContainer, \
    TranslationFunc, \
    FieldInfo, \
    ValidateParameter, \
    ValidateDataField, \
    ValidateProperty, \
    ValidateSpatial, \
    ValidateGridType
from yt.data_objects.yt_array import \
    YTArray
import yt.fields.universal_fields
from yt.fields.particle_fields import \
    particle_deposition_functions, \
    particle_vector_functions, \
    standard_particle_fields
from yt.utilities.physical_constants import \
    mh, \
    mass_sun_cgs
from yt.funcs import *

import yt.utilities.lib as amr_utils

EnzoFieldInfo = FieldInfoContainer.create_with_fallback(FieldInfo, "EFI")
add_field = EnzoFieldInfo.add_field

known_species_masses = dict(
  (sp, mh * v) for sp, v in [
                ("HI": 1.0),
                ("HII": 1.0),
                ("Electron": 1.0),
                ("HeI": 4.0),
                ("HeII": 4.0),
                ("HeIII": 4.0),
                ("H2I": 2.0),
                ("H2II": 2.0),
                ("HM": 1.0),
                ("DI": 2.0),
                ("DII": 2.0),
                ("HDI": 3.0),
    ])

def add_species_field(species, registry):
    # This is currently specific to Enzo.  Hopefully in the future we will have
    # deeper integration with other systems, such as Dengo, to provide better
    # understanding of ionization and molecular states.
    #
    # We have several fields to add based on a given species field.  First off,
    # we add the species field itself.  Then we'll add a few more items...
    #
    registry.add_output_field(("gas", "%s_Density" % species),
                       take_log=True,
                       units="code_mass/code_length**3")
    registry.alias(("gas", "%s_density" % species),
                   ("gas", "%s_Density" % species))
    def _species_mass(field, data):
        return data["gas", "%s_Density"] * data["CellVolume"]
    registry.add_field(("gas", "%s_mass" % species),
                       function=_species_mass,
                       units = "g")
    def _species_fraction(field, data):
        return data["gas", "%s_density"] / data["gas","density"]
    registry.add_field(("gas", "%s_fraction" % species),
                       function=_species_fraction,
                       units = "")
    def _species_number_density(field, data):
        return data["gas", "%s_density"] / _species_mass[species]
    registry.add_field(("gas", "%s_number_density" % species),
                       function=_species_number_density,
                       units = "1/cm**3")

def setup_species_fields(pf, registry, field_list):
    species_names = [field.rsplit("_Density")[0] for field in field_list
                     if field.endswith("_Density")]
    for sp in species_names:
        add_species_field(sp, registry)
    def _number_density(_sp_list, masses):
        def _num_dens_func(field, data):
            num = YTArray(np.zeros(data["Density"].shape, "float64"))
            for sp in _sp_list:
                num += data["%s_density" % sp] / masses[sp]
        return _num_dens_func
    func = _number_density(species_names, known_species_masses)
    registry.add_field(("gas", "number_density"),
                       function = func,
                       units = "1 / cm**3")


b_units = "code_magnetic"
ra_units = "code_length / code_time**2"
rho_units = "code_mass / code_length**3"
vel_units = "code_length / code_time"
_known_other_fields = {
    "Cooling_Time": ("code_time", ["cooling_time"]),
    "HI_kph": ("1/code_time", []),
    "HeI_kph": ("1/code_time", []),
    "HeII_kph": ("1/code_time", []),
    "H2I_kdiss": ("1/code_time", []),
    "Bx": (b_units, ["magnetic_field_x"],
    "By": (b_units, ["magnetic_field_y"],
    "Bz": (b_units, ["magnetic_field_z"],
    "RadAccel1": (ra_units, ["radiation_acceleration_x"],
    "RadAccel2": (ra_units, ["radiation_acceleration_y"],
    "RadAccel3": (ra_units, ["radiation_acceleration_z"],
    "Dark_Matter_Mass": (rho_units, ["dark_matter_mass"]),
    "Temperature": ("K", ["temperature"]),
    "Dust_Temperature": ("K", ["dust_temperature"]),
    "x-velocity": (vel_units, ["velocity_x"]),
    "y-velocity": (vel_units, ["velocity_y"]),
    "z-velocity": (vel_units, ["velocity_z"]),
    "RaySegements": ("", ["ray_segments"]),
    "PhotoGamma": (ra_units, ["photo_gamma"]),
}

_known_particle_fields = {
    "particle_position_x": ("code_length", []),
    "particle_position_y": ("code_length", []),
    "particle_position_z": ("code_length", []),
    "particle_velocity_x": ("code_length / code_time", []),
    "particle_velocity_y": ("code_length / code_time", []),
    "particle_velocity_z": ("code_length / code_time", []),
    "creation_time": ("code_time", []),
    "dynamical_time": ("code_time", []),
    "metallicity_fraction": ("Zsun", []),
    "particle_type": ("", []),
    "particle_index": ("", []),
    "particle_mass": ("code_mass", []),
}

def setup_enzo_gas_fields(registry, field_list, pf):
    # We have a few defaults to set up and alias
    registry.add_output_field(("gas", "Density"),
        units = "code_mass / code_length**3")
    registry.alias(("gas", "density"), ("gas", "Density"))
    registry.add_output_field(("gas", "x-velocity"),
        units = "code_length / code_time")
    registry.add_output_field(("gas", "y-velocity"),
        units = "code_length / code_time")
    registry.add_output_field(("gas", "z-velocity"),
        units = "code_length / code_time")

    # Now we conditionally load a few other things.
    if pf.parameters["MultiSpecies"] > 0:
        setup_species_fields(pf, registry, field_list)
    
    for field in field_list:
        if isinstance(field, tuple): field = field[0]
        args = _known_other_fields.get(field, None)
        if args is None: continue
        units, aliases = args
        registry.add_output_field(
            ("gas", field), units = units)
        for alias in aliases:
            registry.alias(alias, field)

def setup_energy_field(pf, registry, field_list):
    # We check which type of field we need, and then we add it.
    ge_name = None
    te_name = None
    if "Gas_Energy" in field_list:
        ge_name = "Gas_Energy"
    elif "GasEnergy" in field_list:
        ge_name = "GasEnergy"
    if "Total_Energy" in field_list:
        te_name = "Total_Energy"
    elif "TotalEnergy" in field_list:
        te_name = "TotalEnergy"

    if pf.parameters["HydroMethod"] == 2:
        registry.add_output_field(
            ("gas", te_name),
            units="code_length**2/code_time**2"))
        registry.alias(
            ("gas", "thermal_energy"),
            ("gas", te_name))

    elif pf.parameters["DualEnergyFormalism"] == 1:
        registry.add_output_field(
            ("gas", ge_name),
            units="code_length**2/code_time**2")
        registry.alias(
            ("gas", "thermal_energy"),
            ("gas", ge_name),
            units = "erg/g")
    elif pf.paramters["HydroMethod"] in (4, 6):
        registry.add_output_field(
            ("gas", te_name),
            units="code_length**2/code_time**2")
        # Subtract off B-field energy
        def _sub_b(field, data):
            return data[te_name] - 0.5*(
                data["x-velocity"]**2.0
                + data["y-velocity"]**2.0
                + data["z-velocity"]**2.0 ) \
                - data["MagneticEnergy"]/data["Density"]
        registry.add_field(
            ("gas", "thermal_energy"),
            function=_sub_b, units = "erg/g")
    else: # Otherwise, we assume TotalEnergy is kinetic+thermal
        registry.add_output_field(
            ("gas", te_name),
            units = "code_length**2/code_time**2")
        def _tot_minus_kin(field, data):
            return data[te_name] - 0.5*(
                data["x-velocity"]**2.0
                + data["y-velocity"]**2.0
                + data["z-velocity"]**2.0 )
        registry.add_field(
            ("gas", "thermal_energy"),
            units = "erg/g")

    def _kin_energy(field, data):
        return 0.5*data["density"]  * ( data["x-velocity"]**2
                                      + data["y-velocity"]**2.0
                                      + data["z-velocity"]**2.0 )
    registry.add_field("kinetic_energy", function = _kin_energy,
        units = "erg / cm**3")

# Particle functions

# We have now multiplied by grid.dds.prod() inside the IO function.
# So here we multiply just by the conversion to density.

def setup_enzo_particle_fields(registry, ptype):
    for f, (units, aliases) in sorted(_known_particle_fields.items()):
        registry.add_output_field((ptype, f),
            units = units, particle_type = True)
        for alias in aliases:
            registry.alias((alias, (ptype, f))

    particle_vector_functions(ptype,
            ["particle_position_%s" % ax for ax in 'xyz'],
            ["particle_velocity_%s" % ax for ax in 'xyz'],
            registry)
    particle_deposition_functions(ptype, "Coordinates",
        "particle_mass", registry)

    def _age(field, data):
        return data.pf.current_time - data["creation_time"]
    registry.add_field((ptype, "age"), function = _age,
                       particle_type = True,
                       units = "yr")
    standard_particle_fields(registry, ptype)
