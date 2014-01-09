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
from yt.units.yt_array import \
    YTArray

from yt.utilities.physical_constants import \
    mh, \
    mass_sun_cgs

import yt.utilities.lib as amr_utils

b_units = "code_magnetic"
ra_units = "code_length / code_time**2"
rho_units = "code_mass / code_length**3"
vel_units = "code_velocity"

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


class EnzoFieldInfo(FieldInfoContainer):
    known_other_fields = (
        ("Cooling_Time", ("code_time", ["cooling_time"], None)),
        ("HI_kph", ("1/code_time", [], None)),
        ("HeI_kph", ("1/code_time", [], None)),
        ("HeII_kph", ("1/code_time", [], None)),
        ("H2I_kdiss", ("1/code_time", [], None)),
        ("Bx", (b_units, ["magnetic_field_x"], None)),
        ("By", (b_units, ["magnetic_field_y"], None)),
        ("Bz", (b_units, ["magnetic_field_z"], None)),
        ("RadAccel1", (ra_units, ["radiation_acceleration_x"], None)),
        ("RadAccel2", (ra_units, ["radiation_acceleration_y"], None)),
        ("RadAccel3", (ra_units, ["radiation_acceleration_z"], None)),
        ("Dark_Matter_Mass", (rho_units, ["dark_matter_mass"], None)),
        ("Temperature", ("K", ["temperature"], None)),
        ("Dust_Temperature", ("K", ["dust_temperature"], None)),
        ("x-velocity", (vel_units, ["velocity_x"], None)),
        ("y-velocity", (vel_units, ["velocity_y"], None)),
        ("z-velocity", (vel_units, ["velocity_z"], None)),
        ("RaySegments", ("", ["ray_segments"], None)),
        ("PhotoGamma", (ra_units, ["photo_gamma"], None)),
        ("Density", (rho_units, ["density"], None)),
        ("Metal_Density", (rho_units, ["metal_density"], None)),
    )

    known_particle_fields = (
        ("particle_position_x", ("code_length", [], None)),
        ("particle_position_y", ("code_length", [], None)),
        ("particle_position_z", ("code_length", [], None)),
        ("particle_velocity_x", (vel_units, [], None)),
        ("particle_velocity_y", (vel_units, [], None)),
        ("particle_velocity_z", (vel_units, [], None)),
        ("creation_time", ("code_time", [], None)),
        ("dynamical_time", ("code_time", [], None)),
        ("metallicity_fraction", ("code_metallicity", [], None)),
        ("metallicity", ("", [], None)),
        ("particle_type", ("", [], None)),
        ("particle_index", ("", [], None)),
        ("particle_mass", ("code_mass", [], None)),
        ("GridID", ("", [], None)),
        ("identifier", ("", ["particle_index"], None)),
        ("level", ("", [], None)),
    )

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
        elif self.pf.parameters["HydroMethod"] in (4, 6):
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

    def setup_particle_fields(self, ptype):

        def _age(field, data):
            return data.pf.current_time - data["creation_time"]
        self.add_field((ptype, "age"), function = _age,
                           particle_type = True,
                           units = "yr")

        super(EnzoFieldInfo, self).setup_particle_fields(ptype)
