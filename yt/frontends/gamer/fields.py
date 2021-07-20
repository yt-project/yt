import numpy as np

from yt.fields.field_info_container import FieldInfoContainer

b_units = "code_magnetic"
pre_units = "code_mass / (code_length*code_time**2)"
erg_units = "code_mass / (code_length*code_time**2)"
rho_units = "code_mass / code_length**3"
mom_units = "code_mass / (code_length**2*code_time)"
vel_units = "code_velocity"
pot_units = "code_length**2/code_time**2"

psi_units = "code_mass**(1/2) / code_length**(3/2)"


class GAMERFieldInfo(FieldInfoContainer):
    known_other_fields = (
        # hydro fields on disk (GAMER outputs conservative variables)
        ("Dens", (rho_units, [], None)),
        ("MomX", (mom_units, ["momentum_density_x"], None)),
        ("MomY", (mom_units, ["momentum_density_y"], None)),
        ("MomZ", (mom_units, ["momentum_density_z"], None)),
        ("Engy", (erg_units, ["total_energy_density"], None)),
        ("Pote", (pot_units, ["gravitational_potential"], None)),
        # MHD fields on disk (CC=cell-centered)
        ("CCMagX", (b_units, [], "B_x")),
        ("CCMagY", (b_units, [], "B_y")),
        ("CCMagZ", (b_units, [], "B_z")),
        # psiDM fields on disk
        ("Real", (psi_units, ["psidm_real_part"], None)),
        ("Imag", (psi_units, ["psidm_imaginary_part"], None)),
        # particle fields on disk (deposited onto grids)
        ("ParDens", (rho_units, ["particle_density_on_grid"], None)),
        ("TotalDens", (rho_units, ["total_density_on_grid"], None)),
    )

    known_particle_fields = (
        ("ParMass", ("code_mass", ["particle_mass"], None)),
        ("ParPosX", ("code_length", ["particle_position_x"], None)),
        ("ParPosY", ("code_length", ["particle_position_y"], None)),
        ("ParPosZ", ("code_length", ["particle_position_z"], None)),
        ("ParVelX", ("code_velocity", ["particle_velocity_x"], None)),
        ("ParVelY", ("code_velocity", ["particle_velocity_y"], None)),
        ("ParVelZ", ("code_velocity", ["particle_velocity_z"], None)),
        ("ParCreTime", ("code_time", ["particle_creation_time"], None)),
    )

    def __init__(self, ds, field_list):
        super().__init__(ds, field_list)

    # add primitive and other derived variables
    def setup_fluid_fields(self):
        pc = self.ds.units.physical_constants
        from yt.fields.magnetic_field import setup_magnetic_field_aliases

        unit_system = self.ds.unit_system

        if self.ds.srhd:

            c2 = pc.clight * pc.clight

            if self.ds.eos == 4:

                # adiabatic (effective) gamma
                def _gamma(field, data):
                    kT = data["gamer", "Temp"]
                    x = 2.25 * kT / np.sqrt(2.25 * kT * kT + 1.0)
                    c_p = 2.5 + x
                    c_v = 1.5 + x
                    return c_p / c_v

                def htilde(data):
                    kT = data["gamer", "Temp"]
                    x = 2.25 * kT * kT
                    ht = 2.5 * kT + x / (1.0 + np.sqrt(x + 1.0))
                    return ht * c2

                def _sound_speed(field, data):
                    h = htilde(data) / c2 + 1.0
                    kT = data["gamer", "Temp"]
                    cs2 = kT / (3.0 * h)
                    cs2 *= (5.0 * h - 8.0 * kT) / (h - kT)
                    return pc.clight * np.sqrt(cs2)

            else:

                # adiabatic gamma
                def _gamma(field, data):
                    return self.ds.gamma * data["gas", "ones"]

                def htilde(data):
                    kT = data["gamer", "Temp"]
                    g = data["gas", "gamma"]
                    ht = g * kT / (g - 1.0)
                    return ht * c2

                def _sound_speed(field, data):
                    h = htilde(data) / c2 + 1.0
                    cs2 = data["gas", "gamma"] / h * data["gamer", "Temp"]
                    return pc.clight * np.sqrt(cs2)

            # coordinate frame density
            self.alias(
                ("gas", "frame_density"),
                ("gamer", "Dens"),
                units=unit_system["density"],
            )

            self.add_field(
                ("gas", "gamma"), sampling_type="cell", function=_gamma, units=""
            )

            # 4-velocity spatial components
            def four_velocity_xyz(u):
                def _four_velocity(field, data):
                    ui = data["gas", f"momentum_density_{u}"] * c2
                    ui /= data["gas", "frame_density"] * (htilde(data) + c2)
                    return ui

                return _four_velocity

            for u in "xyz":
                self.add_field(
                    ("gas", f"four_velocity_{u}"),
                    sampling_type="cell",
                    function=four_velocity_xyz(u),
                    units=unit_system["velocity"],
                )

            # lorentz factor
            def _lorentz_factor(field, data):
                u2 = data["gas", "four_velocity_magnitude"] ** 2
                return np.sqrt(1.0 + u2 / c2)

            self.add_field(
                ("gas", "lorentz_factor"),
                sampling_type="cell",
                function=_lorentz_factor,
                units="",
            )

            # velocity
            def velocity_xyz(v):
                def _velocity(field, data):
                    return (
                        data["gas", f"four_velocity_{v}"]
                        / data["gas", "lorentz_factor"]
                    )

                return _velocity

            for v in "xyz":
                self.add_field(
                    ("gas", f"velocity_{v}"),
                    sampling_type="cell",
                    function=velocity_xyz(v),
                    units=unit_system["velocity"],
                )

            # density
            def _density(field, data):
                return data["gas", "frame_density"] / data["gas", "lorentz_factor"]

            self.add_field(
                ("gas", "density"),
                sampling_type="cell",
                function=_density,
                units=unit_system["density"],
            )

            # pressure
            def _pressure(field, data):
                return data["gas", "density"] * c2 * data["gamer", "Temp"]

            # thermal energy per mass (i.e., specific)
            def _specific_thermal_energy(field, data):
                eps = data["gas", "density"] * htilde(data) - data["gas", "pressure"]
                return eps / data["gas", "density"]

            # total energy per mass
            def _specific_total_energy(field, data):
                E = data["gamer", "Engy"] + data["gamer", "Dens"] * c2
                return E / data["gamer", "Dens"]

            def _kinetic_energy_density(field, data):
                u2 = data["gas", "four_velocity_magnitude"] ** 2
                gm1 = u2 / c2 / (data["gas", "lorentz_factor"] + 1.0)
                h = htilde(data) + c2
                return gm1 * (data["gamer", "Dens"] * h + data["gas", "pressure"])

            self.add_field(
                ("gas", "kinetic_energy_density"),
                sampling_type="cell",
                function=_kinetic_energy_density,
                units=unit_system["pressure"],
            )

            self.add_field(
                ("gas", "sound_speed"),
                sampling_type="cell",
                function=_sound_speed,
                units=unit_system["velocity"],
            )

            def _mach_number(field, data):
                c_s = data["gas", "sound_speed"]
                u_s = c_s / np.sqrt(1.0 - c_s * c_s / c2)
                return data["gas", "four_velocity_magnitude"] / u_s

            self.add_field(
                ("gas", "mach_number"),
                sampling_type="cell",
                function=_mach_number,
                units="",
            )

        else:

            # density
            self.alias(
                ("gas", "density"),
                ("gamer", "Dens"),
                units=unit_system["density"],
            )

            # velocity
            def velocity_xyz(v):
                def _velocity(field, data):
                    return data["gas", f"momentum_density_{v}"] / data["gas", "density"]

                return _velocity

            for v in "xyz":
                self.add_field(
                    ("gas", f"velocity_{v}"),
                    sampling_type="cell",
                    function=velocity_xyz(v),
                    units=unit_system["velocity"],
                )

            # ====================================================
            # note that yt internal fields assume
            #    [specific_thermal_energy]   = [energy per mass]
            #    [kinetic_energy_density]    = [energy per volume]
            #    [magnetic_energy_density]   = [energy per volume]
            # and we further adopt
            #    [specific_total_energy]     = [energy per mass]
            #    [total_energy_density]      = [energy per volume]
            # ====================================================

            # thermal energy per volume
            def et(data):
                ek = (
                    0.5
                    * (
                        data["gamer", "MomX"] ** 2
                        + data["gamer", "MomY"] ** 2
                        + data["gamer", "MomZ"] ** 2
                    )
                    / data["gamer", "Dens"]
                )
                Et = data["gamer", "Engy"] - ek
                if self.ds.mhd:
                    # magnetic_energy is a yt internal field
                    Et -= data["gas", "magnetic_energy_density"]
                return Et

            # thermal energy per mass (i.e., specific)
            def _specific_thermal_energy(field, data):
                return et(data) / data["gamer", "Dens"]

            # total energy per mass
            def _specific_total_energy(field, data):
                return data["gamer", "Engy"] / data["gamer", "Dens"]

            # pressure
            def _pressure(field, data):
                return et(data) * (data.ds.gamma - 1.0)

        self.add_field(
            ("gas", "specific_thermal_energy"),
            sampling_type="cell",
            function=_specific_thermal_energy,
            units=unit_system["specific_energy"],
        )

        self.add_field(
            ("gas", "specific_total_energy"),
            sampling_type="cell",
            function=_specific_total_energy,
            units=unit_system["specific_energy"],
        )

        self.add_field(
            ("gas", "pressure"),
            sampling_type="cell",
            function=_pressure,
            units=unit_system["pressure"],
        )

        # mean molecular weight
        if hasattr(self.ds, "mu"):

            def _mu(field, data):
                return data.ds.mu * data["index", "ones"]

            self.add_field(
                ("gas", "mean_molecular_weight"),
                sampling_type="cell",
                function=_mu,
                units="",
            )

        # temperature
        def _temperature(field, data):
            return (
                data.ds.mu
                * data["gas", "pressure"]
                * pc.mh
                / (data["gas", "density"] * pc.kb)
            )

        self.add_field(
            ("gas", "temperature"),
            sampling_type="cell",
            function=_temperature,
            units=unit_system["temperature"],
        )

        # magnetic field aliases --> magnetic_field_x/y/z
        if self.ds.mhd:
            setup_magnetic_field_aliases(self, "gamer", [f"CCMag{v}" for v in "XYZ"])

    def setup_particle_fields(self, ptype):
        super().setup_particle_fields(ptype)
