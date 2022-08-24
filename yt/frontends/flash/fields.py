from yt._typing import KnownFieldsT
from yt.fields.field_info_container import FieldInfoContainer

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
en_units = "code_mass * (code_length/code_time)**2"
rho_units = "code_mass / code_length**3"


class FLASHFieldInfo(FieldInfoContainer):
    known_other_fields: KnownFieldsT = (
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
        ("eion", (en_units, [], "Ion Internal Specific Energy")),
        ("eele", (en_units, [], "Electron Internal Specific Energy")),
        ("erad", (en_units, [], "Radiation Internal Specific Energy")),
        ("pden", (rho_units, [], "Particle Mass Density")),
        ("depo", ("code_length**2/code_time**2", [], None)),
        ("ye", ("", [], "Y_e")),
        ("magp", (pres_units, [], None)),
        ("divb", ("code_magnetic/code_length", [], None)),
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
        (
            "particle_gpot",
            ("code_length**2/code_time**2", ["particle_gravitational_potential"], None),
        ),
    )

    def setup_fluid_fields(self):
        from yt.fields.magnetic_field import setup_magnetic_field_aliases

        unit_system = self.ds.unit_system
        # Adopt FLASH 4.6 value for Na
        Na = self.ds.quan(6.022140857e23, "g**-1")
        for i in range(1, 1000):
            self.add_output_field(
                ("flash", f"r{i:03}"),
                sampling_type="cell",
                units="",
                display_name=f"Energy Group {i}",
            )

        # Add energy fields
        def ekin(data):
            ek = data["flash", "velx"] ** 2
            if data.ds.dimensionality >= 2:
                ek += data["flash", "vely"] ** 2
            if data.ds.dimensionality == 3:
                ek += data["flash", "velz"] ** 2
            return 0.5 * ek

        if ("flash", "ener") in self.field_list:
            self.add_output_field(
                ("flash", "ener"),
                sampling_type="cell",
                units="code_length**2/code_time**2",
            )
            self.alias(
                ("gas", "specific_total_energy"),
                ("flash", "ener"),
                units=unit_system["specific_energy"],
            )

        else:

            def _ener(field, data):
                ener = data["flash", "eint"] + ekin(data)
                try:
                    ener += data["flash", "magp"] / data["flash", "dens"]
                except Exception:
                    pass
                return ener

            self.add_field(
                ("gas", "specific_total_energy"),
                sampling_type="cell",
                function=_ener,
                units=unit_system["specific_energy"],
            )
        if ("flash", "eint") in self.field_list:
            self.add_output_field(
                ("flash", "eint"),
                sampling_type="cell",
                units="code_length**2/code_time**2",
            )
            self.alias(
                ("gas", "specific_thermal_energy"),
                ("flash", "eint"),
                units=unit_system["specific_energy"],
            )
        else:

            def _eint(field, data):
                eint = data["flash", "ener"] - ekin(data)
                try:
                    eint -= data["flash", "magp"] / data["flash", "dens"]
                except Exception:
                    pass
                return eint

            self.add_field(
                ("gas", "specific_thermal_energy"),
                sampling_type="cell",
                function=_eint,
                units=unit_system["specific_energy"],
            )

        ## Derived FLASH Fields

        if ("flash", "abar") in self.field_list:
            self.alias(("gas", "mean_molecular_weight"), ("flash", "abar"))
        elif ("flash", "sumy") in self.field_list:

            def _abar(field, data):
                return 1.0 / data["flash", "sumy"]

            self.add_field(
                ("gas", "mean_molecular_weight"),
                sampling_type="cell",
                function=_abar,
                units="",
            )
        elif "eos_singlespeciesa" in self.ds.parameters:

            def _abar(field, data):
                return data.ds.parameters["eos_singlespeciesa"] * data["index", "ones"]

            self.add_field(
                ("gas", "mean_molecular_weight"),
                sampling_type="cell",
                function=_abar,
                units="",
            )

        if ("flash", "sumy") in self.field_list:

            def _nele(field, data):
                return data["flash", "dens"] * data["flash", "ye"] * Na

            self.add_field(
                ("gas", "El_number_density"),
                sampling_type="cell",
                function=_nele,
                units=unit_system["number_density"],
            )

            def _nion(field, data):
                return data["flash", "dens"] * data["flash", "sumy"] * Na

            self.add_field(
                ("gas", "ion_number_density"),
                sampling_type="cell",
                function=_nion,
                units=unit_system["number_density"],
            )

            def _number_density(field, data):
                return (
                    data["gas", "El_number_density"] + data["gas", "ion_number_density"]
                )

        else:

            def _number_density(field, data):
                return data["flash", "dens"] * Na / data["gas", "mean_molecular_weight"]

        self.add_field(
            ("gas", "number_density"),
            sampling_type="cell",
            function=_number_density,
            units=unit_system["number_density"],
        )

        setup_magnetic_field_aliases(self, "flash", [f"mag{ax}" for ax in "xyz"])
