from yt._typing import KnownFieldsT
from yt.fields.field_info_container import FieldInfoContainer
from yt.fields.tensor_fields import setup_stress_energy_ideal

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

        if self.ds.srhd:

            c2 = pc.clight * pc.clight
            c = pc.clight.in_units("code_length / code_time")
            if self.ds.eos == 4:
                fgen = SRHDFields(self.ds.eos, 0.0, c.d)
            else:
                fgen = SRHDFields(self.ds.eos, self.ds.gamma, c.d)

            def _sound_speed(field, data):
                out = fgen.sound_speed(data["gamer", "Temp"].d)
                return data.ds.arr(out, "code_velocity").to(unit_system["velocity"])

            def _gamma(field, data):
                out = fgen.gamma_field(data["gamer", "Temp"].d)
                return data.ds.arr(out, "dimensionless")

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
                    out = fgen.four_velocity_xyz(
                        data["gamer", f"Mom{u.upper()}"].d,
                        data["gamer", "Dens"].d,
                        data["gamer", "Temp"].d,
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
            def _lorentz_factor(field, data):
                out = fgen.lorentz_factor(
                    data["gamer", "Dens"].d,
                    data["gamer", "MomX"].d,
                    data["gamer", "MomY"].d,
                    data["gamer", "MomZ"].d,
                    data["gamer", "Temp"].d,
                )
                return data.ds.arr(out, "dimensionless")

            self.add_field(
                ("gas", "lorentz_factor"),
                sampling_type="cell",
                function=_lorentz_factor,
                units="",
            )

            # velocity
            def velocity_xyz(v):
                def _velocity(field, data):
                    out = fgen.velocity_xyz(
                        data["gamer", "Dens"].d,
                        data["gamer", "MomX"].d,
                        data["gamer", "MomY"].d,
                        data["gamer", "MomZ"].d,
                        data["gamer", "Temp"].d,
                        data["gamer", f"Mom{v.upper()}"].d,
                    )
                    return data.ds.arr(out, "code_velocity").to(unit_system["velocity"])

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
                dens = fgen.density(
                    data["gamer", "Dens"].d,
                    data["gamer", "MomX"].d,
                    data["gamer", "MomY"].d,
                    data["gamer", "MomZ"].d,
                    data["gamer", "Temp"].d,
                )
                return data.ds.arr(dens, rho_units).to(unit_system["density"])

            self.add_field(
                ("gas", "density"),
                sampling_type="cell",
                function=_density,
                units=unit_system["density"],
            )

            # pressure
            def _pressure(field, data):
                out = fgen.pressure(
                    data["gamer", "Dens"].d,
                    data["gamer", "MomX"].d,
                    data["gamer", "MomY"].d,
                    data["gamer", "MomZ"].d,
                    data["gamer", "Temp"].d,
                )
                return data.ds.arr(out, pre_units).to(unit_system["pressure"])

            # thermal energy per mass (i.e., specific)
            def _specific_thermal_energy(field, data):
                out = fgen.specific_thermal_energy(
                    data["gamer", "Dens"].d,
                    data["gamer", "MomX"].d,
                    data["gamer", "MomY"].d,
                    data["gamer", "MomZ"].d,
                    data["gamer", "Temp"].d,
                )
                return data.ds.arr(out, "code_length**2 / code_time**2").to(
                    unit_system["specific_energy"]
                )

            # total energy per mass
            def _specific_total_energy(field, data):
                E = data["gamer", "Engy"] + data["gamer", "Dens"] * c2
                return E / data["gamer", "Dens"]

            def _kinetic_energy_density(field, data):
                out = fgen.kinetic_energy_density(
                    data["gamer", "Dens"].d,
                    data["gamer", "MomX"].d,
                    data["gamer", "MomY"].d,
                    data["gamer", "MomZ"].d,
                    data["gamer", "Temp"].d,
                )
                return data.ds.arr(out, erg_units).to(unit_system["pressure"])

            self.add_field(
                ("gas", "kinetic_energy_density"),
                sampling_type="cell",
                function=_kinetic_energy_density,
                units=unit_system["pressure"],
            )

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
                out = fgen.mach_number(
                    data["gamer", "Dens"].d,
                    data["gamer", "MomX"].d,
                    data["gamer", "MomY"].d,
                    data["gamer", "MomZ"].d,
                    data["gamer", "Temp"].d,
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
