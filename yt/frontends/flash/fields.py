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

from yt.fields.field_info_container import \
    FieldInfoContainer
from yt.utilities.physical_constants import \
    Na

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
        ("velx", ("code_length/code_time", ["velocity_x"], None)),
        ("vely", ("code_length/code_time", ["velocity_y"], None)),
        ("velz", ("code_length/code_time", ["velocity_z"], None)),
        ("dens", ("code_mass/code_length**3", ["density"], None)),
        ("temp", ("code_temperature", ["temperature"], None)),
        ("pres", (pres_units, ["pressure"], None)),
        ("gpot", ("code_length**2/code_time**2", ["gravitational_potential"], None)),
        ("gpol", ("code_length**2/code_time**2", [], None)),
        ("tion", ("code_temperature", [], None)),
        ("tele", ("code_temperature", [], None)),
        ("trad", ("code_temperature", [], None)),
        ("pion", (pres_units, [], None)),
        ("pele", (pres_units, [], "Electron Pressure, P_e")),
        ("prad", (pres_units, [], "Radiation Pressure")),
        ("eion", (erg_units, [], "Ion Internal Energy")),
        ("eele", (erg_units, [], "Electron Internal Energy")),
        ("erad", (erg_units, [], "Radiation Internal Energy")),
        ("pden", (rho_units, [], None)),
        ("depo", ("code_length**2/code_time**2", [], None)),
        ("ye", ("", [], "Y_e")),
        ("magp", (pres_units, [], None)),
        ("divb", ("code_magnetic*code_length", [], None)),
        ("game", ("", [], r"\gamma_e\ \rm{(ratio\ of\ specific\ heats)}")),
        ("gamc", ("", [], r"\gamma_c\ \rm{(ratio\ of\ specific\ heats)}")),
        ("flam", ("", [], None)),
        ("absr", ("", [], "Absorption Coefficient")),
        ("emis", ("", [], "Emissivity")),
        ("cond", ("", [], "Conductivity")),
        ("dfcf", ("", [], "Diffusion Equation Scalar")),
        ("fllm", ("", [], "Flux Limit")),
        ("pipe", ("", [], "P_i/P_e")),
        ("tite", ("", [], "T_i/T_e")),
        ("dbgs", ("", [], "Debug for Shocks")),
        ("cham", ("", [], "Chamber Material Fraction")),
        ("targ", ("", [], "Target Material Fraction")),
        ("sumy", ("", [], None)),
        ("mgdc", ("", [], "Emission Minus Absorption Diffusion Terms")),
        ("magx", (b_units, [], "B_x")),
        ("magy", (b_units, [], "B_y")),
        ("magz", (b_units, [], "B_z")),
    )

    known_particle_fields = (
        ("particle_posx", ("code_length", ["particle_position_x"], None)),
        ("particle_posy", ("code_length", ["particle_position_y"], None)),
        ("particle_posz", ("code_length", ["particle_position_z"], None)),
        ("particle_velx", ("code_length/code_time", ["particle_velocity_x"], None)),
        ("particle_vely", ("code_length/code_time", ["particle_velocity_y"], None)),
        ("particle_velz", ("code_length/code_time", ["particle_velocity_z"], None)),
        ("particle_tag", ("", ["particle_index"], None)),
        ("particle_mass", ("code_mass", ["particle_mass"], None)),
    )

    def setup_fluid_fields(self):
        from yt.fields.magnetic_field import \
            setup_magnetic_field_aliases
        unit_system = self.ds.unit_system
        for i in range(1, 1000):
            self.add_output_field(("flash", "r{0:03}".format(i)), 
                units = "",
                display_name="Energy Group {0}".format(i))
        # Add energy fields
        def ekin(data):
            ek = data["flash","velx"]**2
            if data.ds.dimensionality >= 2:
                ek += data["flash","vely"]**2
            if data.ds.dimensionality == 3:
                ek += data["flash","velz"]**2
            return 0.5*ek
        if ("flash","ener") in self.field_list:
            self.add_output_field(("flash","ener"),
                                  units="code_length**2/code_time**2")
            self.alias(("gas","total_energy"),("flash","ener"),
                       units=unit_system["specific_energy"])
        else:
            def _ener(field, data):
                ener = data["flash","eint"]+ekin(data)
                try:
                    ener += data["flash","magp"]/data["flash","dens"]
                except:
                    pass
                return ener
            self.add_field(("gas","total_energy"), function=_ener,
                           units=unit_system["specific_energy"])
        if ("flash","eint") in self.field_list:
            self.add_output_field(("flash","eint"),
                                  units="code_length**2/code_time**2")
            self.alias(("gas","thermal_energy"),("flash","eint"),
                       units=unit_system["specific_energy"])
        else:
            def _eint(field, data):
                eint = data["flash","ener"]-ekin(data)
                try:
                    eint -= data["flash","magp"]/data["flash","dens"]
                except:
                    pass
                return eint
            self.add_field(("gas","thermal_energy"), function=_eint,
                           units=unit_system["specific_energy"])
        ## Derived FLASH Fields
        def _nele(field, data):
            Na_code = data.ds.quan(Na, '1/code_mass')
            return data["flash","dens"]*data["flash","ye"]*Na_code
        self.add_field(('flash','nele'), function=_nele, units="code_length**-3")
        self.add_field(('flash','edens'), function=_nele, units="code_length**-3")
        def _nion(field, data):
            Na_code = data.ds.quan(Na, '1/code_mass')
            return data["flash","dens"]*data["flash","sumy"]*Na_code
        self.add_field(('flash','nion'), function=_nion, units="code_length**-3")
        
        if ("flash", "abar") in self.field_list:
            self.add_output_field(("flash", "abar"), units="1")
        else:
            def _abar(field, data):
                return 1.0 / data["flash","sumy"]
            self.add_field(("flash","abar"), function=_abar, units="1")

        def _number_density(fields,data):
            return (data["nele"]+data["nion"])
        self.add_field(("gas","number_density"), function=_number_density,
                       units=unit_system["number_density"])

        setup_magnetic_field_aliases(self, "flash", ["mag%s" % ax for ax in "xyz"])


