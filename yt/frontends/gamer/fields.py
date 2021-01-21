from yt.fields.field_info_container import FieldInfoContainer
from yt.utilities.physical_constants import kb, mh, clight

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
        from yt.fields.magnetic_field import setup_magnetic_field_aliases

        unit_system = self.ds.unit_system

<<<<<<< HEAD
        # velocity
        def velocity_xyz(v):
            def _velocity(field, data):
                return data["gas", f"momentum_density_{v}"] / data["gas", "density"]
=======
        if self.ds.special_relativity:
>>>>>>> First stab at some relativistic fields

            c2 = clight*clight

            if self.ds.tm_eos:

                # specific enthalpy
                def _specific_enthalpy(field, data):
                    kT = data["gamer", "Temp"]
                    h = 2.5*kT+np.sqrt(2.25*kT**2+1.0)
                    return h*c2

                # adiabatic (effective) gamma
                def _gamma(field, data):
                    y = data["gas","specific_enthalpy"]/(clight*clight)-1.0
                    x = data["gamer","Temp"]
                    return y/(y-x)

            else:

                # adiabatic gamma
                def _gamma(field, data):
                    return self.ds.gamma.data["gas","ones"]

                # specific enthalpy
                def _specific_enthalpy(field, data):
                    kT = data["gamer", "Temp"]
                    g = data["gas", "gamma"]
                    h = 1.0+g*kT/(g-1.0)
                    return h*c2

            self.add_field(
                ("gas", "specific_enthalpy"),
                sampling_type="cell",
                function=_specific_enthalpy,
                units=unit_system["specific_energy"],
            )

<<<<<<< HEAD
        # ============================================================================
        # note that yt internal fields assume
        #    [specific_thermal_energy]          = [energy per mass]
        #    [kinetic_energy_density]           = [energy per volume]
        #    [magnetic_energy_density]          = [energy per volume]
        # and we further adopt
        #    [specific_total_energy]            = [energy per mass]
        #    [total_energy_density]             = [energy per volume]
        # ============================================================================

        # kinetic energy per volume
        def ek(data):
            return (
                0.5
                * (
                    data["gamer", "MomX"] ** 2
                    + data["gamer", "MomY"] ** 2
                    + data["gamer", "MomZ"] ** 2
=======
            self.add_field(
                ("gas", "gamma"),
                sampling_type="cell",
                function=_gamma,
                units=""
            )

            def four_velocity_xyz(v):
                def _four_velocity(field, data):
                    u = data["gas", f"momentum_{v}"]*2
                    u /= data["gamer", "Dens"]*(c2+data["gas", "specific_enthalpy"])
                    return u
                return _four_velocity

            # 4-velocity spatial components
            for u in "xyz":
                self.add_field(
                    ("gas", f"four_velocity_{u}"),
                    sampling_type="cell",
                    function=four_velocity_xyz(u),
                    units=unit_system["velocity"],
>>>>>>> First stab at some relativistic fields
                )

            # lorentz factor
            def _lorentz_factor(field, data):
                u2 = (data["gas", "four_velocity_x"]**2 +
                      data["gas", "four_velocity_y"]**2 +
                      data["gas", "four_velocity_z"]**2)/c2
                return np.sqrt(1.0+u2)

            self.add_field(
                ("gas", "lorentz_factor"),
                sampling_type="cell",
                function=_lorentz_factor,
                units=""
            )

<<<<<<< HEAD
        # thermal energy per volume
        def et(data):
            Et = data["gamer", "Engy"] - ek(data)
            if self.ds.mhd:
                # magnetic_energy is a yt internal field
                Et -= data["gas", "magnetic_energy_density"]
            return Et

        # thermal energy per mass (i.e., specific)
        def _specific_thermal_energy(field, data):
            return et(data) / data["gamer", "Dens"]

        self.add_field(
            ("gas", "specific_thermal_energy"),
            sampling_type="cell",
            function=_specific_thermal_energy,
            units=unit_system["specific_energy"],
        )

        # total energy per mass
        def _specific_total_energy(field, data):
            return data["gamer", "Engy"] / data["gamer", "Dens"]

        self.add_field(
            ("gas", "specific_total_energy"),
            sampling_type="cell",
            function=_specific_total_energy,
            units=unit_system["specific_energy"],
        )
=======
            # 4-velocity t-component

            def _four_velocity_t(field, data):
                return data["gas", "lorentz_factor"]*clight

            self.add_field(
                ("gas", "four_velocity_t"),
                sampling_type="cell",
                function=_four_velocity_t,
                units=""
            )

            # velocity
            def velocity_xyz(v):
                def _velocity(field, data):
                    return data["gas", f"four_velocity_{v}"] / data["gas", "lorentz_factor"]
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
                return data["gamer", "Dens"] / data["gas", "lorentz_factor"]

            self.add_field(
                ("gas", "density"),
                sampling_type="cell",
                function=_density,
                units=unit_system["density"]
            )

            # pressure
            def _pressure(field, data):
                return data["gas", "density"]*c2*data["gamer", "Temp"]

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
                    return data["gas", f"momentum_{v}"] / data["gas", "density"]
                return _velocity

            for v in "xyz":
                self.add_field(
                    ("gas", f"velocity_{v}"),
                    sampling_type="cell",
                    function=velocity_xyz(v),
                    units=unit_system["velocity"],
                )

            # ============================================================================
            # note that yt internal fields assume
            #    [thermal_energy]          = [energy per mass]
            #    [kinetic_energy]          = [energy per volume]
            #    [magnetic_energy]         = [energy per volume]
            # and we further adopt
            #    [total_energy]            = [energy per mass]
            #    [total_energy_per_volume] = [energy per volume]
            # ============================================================================

            # kinetic energy per volume
            def ek(data):
                return (
                    0.5
                    * (
                        data["gamer", "MomX"] ** 2
                        + data["gamer", "MomY"] ** 2
                        + data["gamer", "MomZ"] ** 2
                    )
                    / data["gamer", "Dens"]
                )

            # thermal energy per volume
            def et(data):
                Et = data["gamer", "Engy"] - ek(data)
                if self.ds.mhd:
                    # magnetic_energy is a yt internal field
                    Et -= data["gas", "magnetic_energy"]
                return Et

            # thermal energy per mass (i.e., specific)
            def _thermal_energy(field, data):
                return et(data) / data["gamer", "Dens"]

            self.add_field(
                ("gas", "thermal_energy"),
                sampling_type="cell",
                function=_thermal_energy,
                units=unit_system["specific_energy"],
            )

            # total energy per mass
            def _total_energy(field, data):
                return data["gamer", "Engy"] / data["gamer", "Dens"]

            self.add_field(
                ("gas", "total_energy"),
                sampling_type="cell",
                function=_total_energy,
                units=unit_system["specific_energy"],
            )
>>>>>>> First stab at some relativistic fields

            # pressure
            def _pressure(field, data):
                return et(data) * (data.ds.gamma - 1.0)

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
                * mh
                / (data["gas", "density"] * kb)
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
