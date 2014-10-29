"""
FLASH-specific fields



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
from yt.fields.field_info_container import \
    FieldInfoContainer
from yt.utilities.physical_constants import \
    kboltz, mh, Na
from yt.units.yt_array import \
    YTArray

# Common fields in FLASH: (Thanks to John ZuHone for this list)
#
# dens gas mass density (g/cc) --
# eint internal energy (ergs/g) --
# ener total energy (ergs/g), with 0.5*v^2 --
# gamc gamma defined as ratio of specific heats, no units
# game gamma defined as in , no units
# gpol gravitational potential from the last timestep (ergs/g)
# gpot gravitational potential from the current timestep (ergs/g)
# grac gravitational acceleration from the current timestep (cm s^-2)
# pden particle mass density (usually dark matter) (g/cc)
# pres pressure (erg/cc)
# temp temperature (K) --
# velx velocity x (cm/s) --
# vely velocity y (cm/s) --
# velz velocity z (cm/s) --

b_units = "code_magnetic"
pres_units = "code_mass/(code_length*code_time**2)"
erg_units = "code_mass * (code_length/code_time)**2"
rho_units = "code_mass / code_length**3"

class FLASHFieldInfo(FieldInfoContainer):
    known_other_fields = (
        (b"velx", ("code_length/code_time", ["velocity_x"], None)),
        (b"vely", ("code_length/code_time", ["velocity_y"], None)),
        (b"velz", ("code_length/code_time", ["velocity_z"], None)),
        (b"dens", ("code_mass/code_length**3", ["density"], None)),
        (b"temp", ("code_temperature", ["temperature"], None)),
        (b"pres", (pres_units, ["pressure"], None)),
        (b"gpot", ("code_length**2/code_time**2", ["gravitational_potential"], None)),
        (b"gpol", ("code_length**2/code_time**2", [], None)),
        (b"tion", ("code_temperature", [], None)),
        (b"tele", ("code_temperature", [], None)),
        (b"trad", ("code_temperature", [], None)),
        (b"pion", (pres_units, [], None)),
        (b"pele", (pres_units, [], "Electron Pressure, P_e")),
        (b"prad", (pres_units, [], "Radiation Pressure")),
        (b"eion", (erg_units, [], "Ion Internal Energy")),
        (b"eele", (erg_units, [], "Electron Internal Energy")),
        (b"erad", (erg_units, [], "Radiation Internal Energy")),
        (b"pden", (rho_units, [], None)),
        (b"depo", ("code_length**2/code_time**2", [], None)),
        (b"ye", ("", [], "Y_e")),
        (b"magp", (pres_units, [], None)),
        (b"divb", ("code_magnetic*code_length", [], None)),
        (b"game", ("", [], "\gamma_e\/\rm{(ratio\/of\/specific\/heats)}")),
        (b"gamc", ("", [], "\gamma_c\/\rm{(ratio\/of\/specific\/heats)}")),
        (b"flam", ("", [], None)),
        (b"absr", ("", [], "Absorption Coefficient")),
        (b"emis", ("", [], "Emissivity")),
        (b"cond", ("", [], "Conductivity")),
        (b"dfcf", ("", [], "Diffusion Equation Scalar")),
        (b"fllm", ("", [], "Flux Limit")),
        (b"pipe", ("", [], "P_i/P_e")),
        (b"tite", ("", [], "T_i/T_e")),
        (b"dbgs", ("", [], "Debug for Shocks")),
        (b"cham", ("", [], "Chamber Material Fraction")),
        (b"targ", ("", [], "Target Material Fraction")),
        (b"sumy", ("", [], None)),
        (b"mgdc", ("", [], "Emission Minus Absorption Diffusion Terms")),
        (b"magx", (b_units, ["magnetic_field_x"], "B_x")),
        (b"magy", (b_units, ["magnetic_field_y"], "B_y")),
        (b"magz", (b_units, ["magnetic_field_z"], "B_z")),
    )

    known_particle_fields = (
        (b"particle_posx", ("code_length", ["particle_position_x"], None)),
        (b"particle_posy", ("code_length", ["particle_position_y"], None)),
        (b"particle_posz", ("code_length", ["particle_position_z"], None)),
        (b"particle_velx", ("code_length/code_time", ["particle_velocity_x"], None)),
        (b"particle_vely", ("code_length/code_time", ["particle_velocity_y"], None)),
        (b"particle_velz", ("code_length/code_time", ["particle_velocity_z"], None)),
        (b"particle_tag", ("", ["particle_index"], None)),
        (b"particle_mass", ("code_mass", ["particle_mass"], None)),
    )

    def setup_fluid_fields(self):
        for i in range(1, 1000):
            self.add_output_field(("flash", "r{0:03}".format(i)), 
                units = "",
                display_name="Energy Group {0}".format(i))
        # Add energy fields
        def ekin(data):
            ek = data["flash",b"velx"]**2
            if data.ds.dimensionality >= 2:
                ek += data["flash",b"vely"]**2
            if data.ds.dimensionality == 3:
                ek += data["flash",b"velz"]**2
            return 0.5*ek
        if ("flash",b"ener") in self.field_list:
            self.add_output_field(("flash",b"ener"),
                                  units="code_length**2/code_time**2")
            self.alias(("gas","total_energy"),("flash",b"ener"),
                       units="erg/g")
        else:
            def _ener(field, data):
                ener = data["flash",b"eint"]+ekin(data)
                try:
                    ener += data["flash",b"magp"]/data["flash",b"dens"]
                except:
                    pass
                return ener
            self.add_field(("gas","total_energy"), function=_ener,
                           units="erg/g")
        if ("flash",b"eint") in self.field_list:
            self.add_output_field(("flash",b"eint"),
                                  units="code_length**2/code_time**2")
            self.alias(("gas","thermal_energy"),("flash",b"eint"),
                       units="erg/g")
        else:
            def _eint(field, data):
                eint = data["flash",b"ener"]-ekin(data)
                try:
                    eint -= data["flash",b"magp"]/data["flash",b"dens"]
                except:
                    pass
                return eint
            self.add_field(("gas","thermal_energy"), function=_eint,
                           units="erg/g")
        ## Derived FLASH Fields
        def _nele(field, data):
            Na_code = data.ds.quan(Na, '1/code_mass')
            return data["flash",b"dens"]*data["flash",b"ye"]*Na_code
        self.add_field(('flash',b'nele'), function=_nele, units="code_length**-3")
        self.add_field(('flash',b'edens'), function=_nele, units="code_length**-3")
        def _nion(field, data):
            Na_code = data.ds.quan(Na, '1/code_mass')
            return data["flash",b"dens"]*data["flash",b"sumy"]*Na_code
        self.add_field(('flash',b'nion'), function=_nion, units="code_length**-3")
        def _abar(field, data):
            try:
                return data["flash",b"abar"]
            except:
                return 1.0/data["flash",b"sumy"]
        self.add_field(("flash",b"abar"), function=_abar, units="1")
        def _number_density(fields,data) :
            return (data[b"nele"]+data[b"nion"])
        self.add_field(("gas","number_density"), function=_number_density,
                       units="cm**-3")


