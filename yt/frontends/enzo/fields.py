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

from yt.funcs import mylog
from yt.fields.field_info_container import \
    FieldInfoContainer
from yt.data_objects.yt_array import \
    YTArray
from yt.fields.particle_fields import \
    particle_deposition_functions, \
    particle_vector_functions, \
    standard_particle_fields
from yt.fields.vector_operations import \
    create_vector_fields

from yt.utilities.physical_constants import \
    mh, \
    mass_sun_cgs

import yt.utilities.lib as amr_utils

b_units = "" # Gotta fix this
ra_units = "code_length / code_time**2"
rho_units = "code_mass / code_length**3"
vel_units = "code_length / code_time"

known_species_masses = dict(
  (sp, mh * v) for sp, v in [
                ("HI", 1.0),
                ("HII", 1.0),
                ("Electron", 1.0),
                ("HeI", 4.0),
                ("HeII", 4.0),
                ("HeIII", 4.0),
                ("H2I", 2.0),
                ("H2II", 2.0),
                ("HM", 1.0),
                ("DI", 2.0),
                ("DII", 2.0),
                ("HDI", 3.0),
    ])

known_other_fields = {
    "Cooling_Time": ("code_time", ["cooling_time"]),
    "HI_kph": ("1/code_time", []),
    "HeI_kph": ("1/code_time", []),
    "HeII_kph": ("1/code_time", []),
    "H2I_kdiss": ("1/code_time", []),
    "Bx": (b_units, ["magnetic_field_x"]),
    "By": (b_units, ["magnetic_field_y"]),
    "Bz": (b_units, ["magnetic_field_z"]),
    "RadAccel1": (ra_units, ["radiation_acceleration_x"]),
    "RadAccel2": (ra_units, ["radiation_acceleration_y"]),
    "RadAccel3": (ra_units, ["radiation_acceleration_z"]),
    "Dark_Matter_Mass": (rho_units, ["dark_matter_mass"]),
    "Temperature": ("K", ["temperature"]),
    "Dust_Temperature": ("K", ["dust_temperature"]),
    "x-velocity": (vel_units, ["velocity_x"]),
    "y-velocity": (vel_units, ["velocity_y"]),
    "z-velocity": (vel_units, ["velocity_z"]),
    "RaySegements": ("", ["ray_segments"]),
    "PhotoGamma": (ra_units, ["photo_gamma"]),
    "Density": (rho_units, ["density"]),
}

known_particle_fields = {
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

class EnzoFieldInfo(FieldInfoContainer):

    def __init__(self, pf, field_list):
        if pf.parameters["HydroMethod"] == 2:
            sl_left = slice(None,-2,None)
            sl_right = slice(1,-1,None)
            div_fac = 1.0
        else:
            sl_left = slice(None,-2,None)
            sl_right = slice(2,None,None)
            div_fac = 2.0
        slice_info = (sl_left, sl_right, div_fac)
        super(EnzoFieldInfo, self).__init__(pf, field_list, slice_info)

    def add_species_field(self, species):
        # This is currently specific to Enzo.  Hopefully in the future we will
        # have deeper integration with other systems, such as Dengo, to provide
        # better understanding of ionization and molecular states.
        #
        # We have several fields to add based on a given species field.  First
        # off, we add the species field itself.  Then we'll add a few more
        # items...
        #
        self.add_output_field(("enzo", "%s_Density" % species),
                           take_log=True,
                           units="code_mass/code_length**3")
        self.alias(("gas", "%s_density" % species),
                   ("enzo", "%s_Density" % species))
        def _species_mass(field, data):
            return data["gas", "%s_density" % species] \
                 * data["cell_volume"]
        self.add_field(("gas", "%s_mass" % species),
                           function=_species_mass,
                           units = "g")
        def _species_fraction(field, data):
            return data["gas", "%s_density" % species] \
                 / data["gas","density"]
        self.add_field(("gas", "%s_fraction" % species),
                           function=_species_fraction,
                           units = "")
        def _species_number_density(field, data):
            return data["gas", "%s_density" % species] \
                / known_species_masses[species]
        self.add_field(("gas", "%s_number_density" % species),
                           function=_species_number_density,
                           units = "1/cm**3")

    def setup_species_fields(self):
        species_names = [fn.rsplit("_Density")[0] for ft, fn in 
                         self.field_list if fn.endswith("_Density")]
        species_names = [sp for sp in species_names
                         if sp in known_species_masses]
        for sp in species_names:
            self.add_species_field(sp)
        def _number_density(_sp_list, masses):
            def _num_dens_func(field, data):
                num = YTArray(np.zeros(data["density"].shape, "float64"))
                for sp in _sp_list:
                    num += data["%s_density" % sp] / masses[sp]
            return _num_dens_func
        func = _number_density(species_names, known_species_masses)
        self.add_field(("gas", "number_density"),
                           function = func,
                           units = "1 / cm**3")

    def setup_fluid_fields(self):
        for field in sorted(self.field_list):
            if not isinstance(field, tuple):
                raise RuntimeError
            args = known_other_fields.get(field[1], None)
            if args is None: continue
            units, aliases = args
            self.add_output_field(field, units = units)
            for alias in aliases:
                self.alias(("gas", alias), field)

        # Now we conditionally load a few other things.
        if self.pf.parameters["MultiSpecies"] > 0:
            self.setup_species_fields()
        self.setup_energy_field()

    def setup_energy_field(self):
        # We check which type of field we need, and then we add it.
        ge_name = None
        te_name = None
        if ("enzo", "Gas_Energy") in self.field_list:
            ge_name = "Gas_Energy"
        elif ("enzo", "GasEnergy") in self.field_list:
            ge_name = "GasEnergy"
        if ("enzo", "Total_Energy") in self.field_list:
            te_name = "Total_Energy"
        elif ("enzo", "TotalEnergy") in self.field_list:
            te_name = "TotalEnergy"

        if self.pf.parameters["HydroMethod"] == 2:
            self.add_output_field(("enzo", te_name),
                units="code_length**2/code_time**2")
            self.alias(("gas", "thermal_energy"), ("enzo", te_name))

        elif self.pf.parameters["DualEnergyFormalism"] == 1:
            self.add_output_field(
                ("enzo", ge_name),
                units="code_length**2/code_time**2")
            self.alias(
                ("gas", "thermal_energy"),
                ("enzo", ge_name),
                units = "erg/g")
        elif self.pf.paramters["HydroMethod"] in (4, 6):
            self.add_output_field(
                ("enzo", te_name),
                units="code_length**2/code_time**2")
            # Subtract off B-field energy
            def _sub_b(field, data):
                return data[te_name] - 0.5*(
                    data["x-velocity"]**2.0
                    + data["y-velocity"]**2.0
                    + data["z-velocity"]**2.0 ) \
                    - data["MagneticEnergy"]/data["Density"]
            self.add_field(
                ("gas", "thermal_energy"),
                function=_sub_b, units = "erg/g")
        else: # Otherwise, we assume TotalEnergy is kinetic+thermal
            self.add_output_field(
                ("enzo", te_name),
                units = "code_length**2/code_time**2")
            def _tot_minus_kin(field, data):
                return data[te_name] - 0.5*(
                    data["x-velocity"]**2.0
                    + data["y-velocity"]**2.0
                    + data["z-velocity"]**2.0 )
            self.add_field(
                ("gas", "thermal_energy"),
                units = "erg/g")

        def _kin_energy(field, data):
            return 0.5*data["density"] * ( data["x-velocity"]**2
                                         + data["y-velocity"]**2.0
                                         + data["z-velocity"]**2.0 )
        self.add_field(("gas", "kinetic_energy"),
            function = _kin_energy,
            units = "erg / cm**3")

    def setup_particle_fields(self, ptype):
        for f, (units, aliases) in sorted(known_particle_fields.items()):
            self.add_output_field((ptype, f),
                units = units, particle_type = True)
            for alias in aliases:
                self.alias(alias, (ptype, f))

        particle_vector_functions(ptype,
                ["particle_position_%s" % ax for ax in 'xyz'],
                ["particle_velocity_%s" % ax for ax in 'xyz'],
                self)
        particle_deposition_functions(ptype, "Coordinates",
            "particle_mass", self)

        def _age(field, data):
            return data.pf.current_time - data["creation_time"]
        self.add_field((ptype, "age"), function = _age,
                           particle_type = True,
                           units = "yr")
        standard_particle_fields(self, ptype)
