from yt._typing import KnownFieldsT
from yt.fields.field_info_container import FieldInfoContainer
from yt.fields.tensor_fields import setup_stress_energy_ideal
from yt.funcs import mylog

from .cfields import SRHDFields

b_units = "code_magnetic"
pre_units = "code_mass / (code_length*code_time**2)"
erg_units = "code_mass / (code_length*code_time**2)"
rho_units = "code_mass / code_length**3"
mom_units = "code_mass / (code_length**2*code_time)"
vel_units = "code_velocity"
pot_units = "code_length**2/code_time**2"

psi_units = "code_mass**(1/2) / code_length**(3/2)"


class GAMERFieldInfo(FieldInfoContainer):
    known_other_fields: KnownFieldsT = (
        # hydro fields on disk (GAMER outputs conservative variables)
        ("Dens", (rho_units, [], None)),
        ("MomX", (mom_units, ["momentum_density_x"], None)),
        ("MomY", (mom_units, ["momentum_density_y"], None)),
        ("MomZ", (mom_units, ["momentum_density_z"], None)),
        ("Engy", (erg_units, [], None)),
        ("CRay", (erg_units, ["cosmic_ray_energy_density"], None)),
        ("Pote", (pot_units, ["gravitational_potential"], None)),
        ("Pres", (pre_units, ["pressure"], None)),
        ("Temp", ("code_temperature", ["temperature"], None)),
        ("Enth", (pot_units, ["specific_reduced_enthalpy"], None)),
        ("Mach", ("dimensionless", ["mach_number"], None)),
        ("Cs", (vel_units, ["sound_speed"], None)),
        ("DivVel", ("1/code_time", ["velocity_divergence"], None)),
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

    known_particle_fields: KnownFieldsT = (
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
        unit_system.registry = self.ds.unit_registry  # TODO: Why do I need this?!

        if self.ds.opt_unit:
            temp_conv = pc.kb / (self.ds.mu * pc.mh)
        else:
            temp_conv = (
                self.ds.arr(1.0, "code_velocity**2/code_temperature") / self.ds.mu
            )

        if self.ds.srhd:
            if self.ds.opt_unit:
                c2 = pc.clight * pc.clight
            else:
                c2 = self.ds.arr(1.0, "code_velocity**2")
            invc2 = 1.0 / c2

            if ("gamer", "Temp") not in self.field_list:
                mylog.warning(
                    'The temperature field "Temp" is not present in the dataset. Most '
                    'SRHD fields will not be available!! Please set "OPT__OUTPUT_TEMP '
                    '= 1" in Input__Parameter and re-run the simulation.'
                )
            if ("gamer", "Enth") not in self.field_list:
                mylog.warning(
                    'The reduced enthalpy field "Enth" is not present in the dataset. '
                    "Most SRHD fields will not be available!! Please set "
                    '"OPT__OUTPUT_ENTHALPY = 1" in Input__Parameter and re-run the '
                    "simulation."
                )

            # EOS functions
            gamma = self.ds.gamma if self.ds.eos == 1 else 0.0
            fgen = SRHDFields(self.ds.eos, gamma)

            # temperature fraction (kT/mc^2)
            def _temp_fraction(field, data):
                return data["gamer", "Temp"] * temp_conv * invc2

            self.add_field(
                ("gas", "temp_fraction"),
                function=_temp_fraction,
                sampling_type="cell",
                units="",
            )

            # specific enthalpy
            def _specific_enthalpy(field, data):
                return data["gas", "specific_reduced_enthalpy"] + c2

            self.add_field(
                ("gas", "specific_enthalpy"),
                function=_specific_enthalpy,
                sampling_type="cell",
                units=unit_system["specific_energy"],
            )

            # sound speed
            if ("gamer", "Cs") not in self.field_list:

                def _sound_speed(field, data):
                    out = fgen.sound_speed(
                        data["gas", "temp_fraction"].d,
                        data["gamer", "Enth"].d,
                    )
                    return data.ds.arr(out, "code_velocity").to(unit_system["velocity"])

                self.add_field(
                    ("gas", "sound_speed"),
                    sampling_type="cell",
                    function=_sound_speed,
                    units=unit_system["velocity"],
                )

            # ratio of specific heats (gamma)
            def _gamma(field, data):
                out = fgen.gamma_field(data["gas", "temp_fraction"].d)
                return data.ds.arr(out, "dimensionless")

            self.add_field(
                ("gas", "gamma"), sampling_type="cell", function=_gamma, units=""
            )

            # reduced total energy density
            self.alias(
                ("gas", "reduced_total_energy_density"),
                ("gamer", "Engy"),
                units=unit_system["pressure"],
            )

            # total energy density
            def _total_energy_density(field, data):
                return data["gamer", "Engy"] + data["gamer", "Dens"] * c2

            self.add_field(
                ("gas", "total_energy_density"),
                sampling_type="cell",
                function=_total_energy_density,
                units=unit_system["pressure"],
            )

            # coordinate frame density
            self.alias(
                ("gas", "frame_density"),
                ("gamer", "Dens"),
                units=unit_system["density"],
            )

            # 4-velocity spatial components
            def four_velocity_xyz(u):
                def _four_velocity(field, data):
                    out = fgen.four_velocity_xyz(
                        data["gamer", "Dens"].d,
                        data["gamer", f"Mom{u.upper()}"].d,
                        data["gamer", "Enth"].d,
                    )
                    return data.ds.arr(out, "code_velocity").to(unit_system["velocity"])

                return _four_velocity

            for u in "xyz":
                self.add_field(
                    ("gas", f"four_velocity_{u}"),
                    sampling_type="cell",
                    function=four_velocity_xyz(u),
                    units=unit_system["velocity"],
                )

            # lorentz factor
            if ("gamer", "Lrtz") in self.field_list:

                def _lorentz_factor(field, data):
                    return data["gamer", "Lrtz"]

            else:

                def _lorentz_factor(field, data):
                    out = fgen.lorentz_factor(
                        data["gamer", "Dens"].d,
                        data["gamer", "MomX"].d,
                        data["gamer", "MomY"].d,
                        data["gamer", "MomZ"].d,
                        data["gamer", "Enth"].d,
                    )
                    return data.ds.arr(out, "dimensionless")

            self.add_field(
                ("gas", "lorentz_factor"),
                sampling_type="cell",
                function=_lorentz_factor,
                units="",
            )

            # density
            def _density(field, data):
                return data["gamer", "Dens"] / data["gas", "lorentz_factor"]

            self.add_field(
                ("gas", "density"),
                sampling_type="cell",
                function=_density,
                units=unit_system["density"],
            )

            # pressure
            def _pressure(field, data):
                p = data["gas", "density"] * data["gas", "temp_fraction"]
                return p * c2

            # thermal energy per mass (i.e., specific)
            def _specific_thermal_energy(field, data):
                return (
                    data["gas", "specific_reduced_enthalpy"]
                    - c2 * data["gas", "temp_fraction"]
                )

            # total energy per mass
            def _specific_total_energy(field, data):
                return data["gas", "total_energy_density"] / data["gas", "density"]

            # kinetic energy density
            def _kinetic_energy_density(field, data):
                out = fgen.kinetic_energy_density(
                    data["gamer", "Dens"].d,
                    data["gamer", "MomX"].d,
                    data["gamer", "MomY"].d,
                    data["gamer", "MomZ"].d,
                    data["gas", "temp_fraction"].d,
                    data["gamer", "Enth"].d,
                )
                return data.ds.arr(out, erg_units).to(unit_system["pressure"])

            self.add_field(
                ("gas", "kinetic_energy_density"),
                sampling_type="cell",
                function=_kinetic_energy_density,
                units=unit_system["pressure"],
            )

            # Mach number
            if ("gamer", "Mach") not in self.field_list:

                def _mach_number(field, data):
                    out = fgen.mach_number(
                        data["gamer", "Dens"].d,
                        data["gamer", "MomX"].d,
                        data["gamer", "MomY"].d,
                        data["gamer", "MomZ"].d,
                        data["gas", "temp_fraction"].d,
                        data["gamer", "Enth"].d,
                    )
                    return data.ds.arr(out, "dimensionless")

                self.add_field(
                    ("gas", "mach_number"),
                    sampling_type="cell",
                    function=_mach_number,
                    units="",
                )

            setup_stress_energy_ideal(self)

        else:  # not RHD
            # density
            self.alias(
                ("gas", "density"), ("gamer", "Dens"), units=unit_system["density"]
            )

            self.alias(
                ("gas", "total_energy_density"),
                ("gamer", "Engy"),
                units=unit_system["pressure"],
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
                if getattr(self.ds, "gamma_cr", None):
                    # cosmic rays are included in this dataset
                    Et -= data["gas", "cosmic_ray_energy_density"]
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

        # velocity
        def velocity_xyz(v):
            if ("gamer", f"Vel{v.upper()}") in self.field_list:

                def _velocity(field, data):
                    return data.ds.arr(
                        data["gamer", f"Vel{v.upper()}"].d, "code_velocity"
                    ).to(unit_system["velocity"])

            elif self.ds.srhd:

                def _velocity(field, data):
                    return (
                        data["gas", f"four_velocity_{v}"]
                        / data["gas", "lorentz_factor"]
                    )

            else:

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

        if ("gamer", "Pres") not in self.field_list:
            self.add_field(
                ("gas", "pressure"),
                sampling_type="cell",
                function=_pressure,
                units=unit_system["pressure"],
            )

        self.add_field(
            ("gas", "specific_thermal_energy"),
            sampling_type="cell",
            function=_specific_thermal_energy,
            units=unit_system["specific_energy"],
        )

        def _thermal_energy_density(field, data):
            return data["gas", "density"] * data["gas", "specific_thermal_energy"]

        self.add_field(
            ("gas", "thermal_energy_density"),
            sampling_type="cell",
            function=_thermal_energy_density,
            units=unit_system["pressure"],
        )

        self.add_field(
            ("gas", "specific_total_energy"),
            sampling_type="cell",
            function=_specific_total_energy,
            units=unit_system["specific_energy"],
        )

        if getattr(self.ds, "gamma_cr", None):

            def _cr_pressure(field, data):
                return (data.ds.gamma_cr - 1.0) * data[
                    "gas", "cosmic_ray_energy_density"
                ]

            self.add_field(
                ("gas", "cosmic_ray_pressure"),
                _cr_pressure,
                sampling_type="cell",
                units=self.ds.unit_system["pressure"],
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
        if ("gamer", "Temp") not in self.field_list:

            def _temperature(field, data):
                return data["gas", "pressure"] / (data["gas", "density"] * temp_conv)

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
