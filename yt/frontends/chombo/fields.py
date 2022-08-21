import numpy as np

from yt._typing import KnownFieldsT
from yt.fields.field_info_container import (
    FieldInfoContainer,
    particle_deposition_functions,
    particle_vector_functions,
    standard_particle_fields,
)
from yt.units.unit_object import Unit  # type: ignore
from yt.utilities.exceptions import YTFieldNotFound

rho_units = "code_mass / code_length**3"
mom_units = "code_mass / (code_time * code_length**2)"
eden_units = "code_mass / (code_time**2 * code_length)"  # erg / cm^3
vel_units = "code_length / code_time"
b_units = "code_magnetic"


class ChomboFieldInfo(FieldInfoContainer):
    # no custom behaviour is needed yet
    pass


# Orion 2 Fields
# We duplicate everything here from Boxlib, because we want to be able to
# subclass it and that can be somewhat tricky.
class Orion2FieldInfo(ChomboFieldInfo):
    known_other_fields: KnownFieldsT = (
        ("density", (rho_units, ["density"], None)),
        ("energy-density", (eden_units, ["total_energy_density"], None)),
        ("radiation-energy-density", (eden_units, ["radiation_energy_density"], None)),
        ("X-momentum", (mom_units, ["momentum_density_x"], None)),
        ("Y-momentum", (mom_units, ["momentum_density_y"], None)),
        ("Z-momentum", (mom_units, ["momentum_density_z"], None)),
        ("temperature", ("K", ["temperature"], None)),
        ("X-magnfield", (b_units, [], None)),
        ("Y-magnfield", (b_units, [], None)),
        ("Z-magnfield", (b_units, [], None)),
        ("directrad-dedt-density", (eden_units, ["directrad-dedt-density"], None)),
        ("directrad-dpxdt-density", (mom_units, ["directrad-dpxdt-density"], None)),
        ("directrad-dpydt-density", (mom_units, ["directrad-dpydt-density"], None)),
        ("directrad-dpzdt-density", (mom_units, ["directrad-dpzdt-density"], None)),
    )
    known_particle_fields: KnownFieldsT = (
        ("particle_mass", ("code_mass", [], None)),
        ("particle_position_x", ("code_length", [], None)),
        ("particle_position_y", ("code_length", [], None)),
        ("particle_position_z", ("code_length", [], None)),
        ("particle_momentum_x", ("code_mass*code_length/code_time", [], None)),
        ("particle_momentum_y", ("code_mass*code_length/code_time", [], None)),
        ("particle_momentum_z", ("code_mass*code_length/code_time", [], None)),
        # Note that these are *internal* agmomen
        ("particle_angmomen_x", ("code_length**2/code_time", [], None)),
        ("particle_angmomen_y", ("code_length**2/code_time", [], None)),
        ("particle_angmomen_z", ("code_length**2/code_time", [], None)),
        ("particle_mlast", ("code_mass", [], None)),
        ("particle_r", ("code_length", [], None)),
        ("particle_mdeut", ("code_mass", [], None)),
        ("particle_n", ("", [], None)),
        ("particle_mdot", ("code_mass/code_time", [], None)),
        ("particle_burnstate", ("", [], None)),
        ("particle_luminosity", ("", [], None)),
        ("particle_id", ("", ["particle_index"], None)),
    )

    def setup_particle_fields(self, ptype):
        def _get_vel(axis):
            def velocity(field, data):
                return (
                    data[(ptype, f"particle_momentum_{axis}")]
                    / data[(ptype, "particle_mass")]
                )

            return velocity

        for ax in "xyz":
            self.add_field(
                (ptype, f"particle_velocity_{ax}"),
                sampling_type="particle",
                function=_get_vel(ax),
                units="code_length/code_time",
            )

        super().setup_particle_fields(ptype)

    def setup_fluid_fields(self):
        from yt.fields.magnetic_field import setup_magnetic_field_aliases

        unit_system = self.ds.unit_system

        def _thermal_energy_density(field, data):
            try:
                return (
                    data[("chombo", "energy-density")]
                    - data[("gas", "kinetic_energy_density")]
                    - data[("gas", "magnetic_energy_density")]
                )
            except YTFieldNotFound:
                return (
                    data[("chombo", "energy-density")]
                    - data[("gas", "kinetic_energy_density")]
                )

        def _specific_thermal_energy(field, data):
            return data[("gas", "thermal_energy_density")] / data[("gas", "density")]

        def _magnetic_energy_density(field, data):
            ret = data[("chombo", "X-magnfield")] ** 2
            if data.ds.dimensionality > 1:
                ret = ret + data[("chombo", "Y-magnfield")] ** 2
            if data.ds.dimensionality > 2:
                ret = ret + data[("chombo", "Z-magnfield")] ** 2
            return ret / 8.0 / np.pi

        def _specific_magnetic_energy(field, data):
            return data[("gas", "specific_magnetic_energy")] / data[("gas", "density")]

        def _kinetic_energy_density(field, data):
            p2 = data[("chombo", "X-momentum")] ** 2
            if data.ds.dimensionality > 1:
                p2 = p2 + data[("chombo", "Y-momentum")] ** 2
            if data.ds.dimensionality > 2:
                p2 = p2 + data[("chombo", "Z-momentum")] ** 2
            return 0.5 * p2 / data[("gas", "density")]

        def _specific_kinetic_energy(field, data):
            return data[("gas", "kinetic_energy_density")] / data[("gas", "density")]

        def _temperature(field, data):
            c_v = data.ds.quan(data.ds.parameters["radiation.const_cv"], "erg/g/K")
            return data[("gas", "specific_thermal_energy")] / c_v

        def _get_vel(axis):
            def velocity(field, data):
                return (
                    data[("gas", f"momentum_density_{axis}")] / data[("gas", "density")]
                )

            return velocity

        for ax in "xyz":
            self.add_field(
                ("gas", f"velocity_{ax}"),
                sampling_type="cell",
                function=_get_vel(ax),
                units=unit_system["velocity"],
            )
        self.add_field(
            ("gas", "specific_thermal_energy"),
            sampling_type="cell",
            function=_specific_thermal_energy,
            units=unit_system["specific_energy"],
        )
        self.add_field(
            ("gas", "thermal_energy_density"),
            sampling_type="cell",
            function=_thermal_energy_density,
            units=unit_system["pressure"],
        )
        self.add_field(
            ("gas", "kinetic_energy_density"),
            sampling_type="cell",
            function=_kinetic_energy_density,
            units=unit_system["pressure"],
        )
        self.add_field(
            ("gas", "specific_kinetic_energy"),
            sampling_type="cell",
            function=_specific_kinetic_energy,
            units=unit_system["specific_energy"],
        )
        self.add_field(
            ("gas", "magnetic_energy_density"),
            sampling_type="cell",
            function=_magnetic_energy_density,
            units=unit_system["pressure"],
        )
        self.add_field(
            ("gas", "specific_magnetic_energy"),
            sampling_type="cell",
            function=_specific_magnetic_energy,
            units=unit_system["specific_energy"],
        )
        self.add_field(
            ("gas", "temperature"),
            sampling_type="cell",
            function=_temperature,
            units=unit_system["temperature"],
        )

        setup_magnetic_field_aliases(
            self, "chombo", [f"{ax}-magnfield" for ax in "XYZ"]
        )


class ChomboPICFieldInfo3D(FieldInfoContainer):
    known_other_fields: KnownFieldsT = (
        ("density", (rho_units, ["density", "Density"], None)),
        (
            "potential",
            ("code_length**2 / code_time**2", ["potential", "Potential"], None),
        ),
        ("gravitational_field_x", ("code_length / code_time**2", [], None)),
        ("gravitational_field_y", ("code_length / code_time**2", [], None)),
        ("gravitational_field_z", ("code_length / code_time**2", [], None)),
    )
    known_particle_fields: KnownFieldsT = (
        ("particle_mass", ("code_mass", [], None)),
        ("particle_position_x", ("code_length", [], None)),
        ("particle_position_y", ("code_length", [], None)),
        ("particle_position_z", ("code_length", [], None)),
        ("particle_velocity_x", ("code_length / code_time", [], None)),
        ("particle_velocity_y", ("code_length / code_time", [], None)),
        ("particle_velocity_z", ("code_length / code_time", [], None)),
    )

    # I am re-implementing this here to override a few default behaviors:
    # I don't want to skip output units for code_length and I do want
    # particle_fields to default to take_log = False.
    def setup_particle_fields(self, ptype, ftype="gas", num_neighbors=64):
        skip_output_units = ()
        for f, (units, aliases, dn) in sorted(self.known_particle_fields):
            units = self.ds.field_units.get((ptype, f), units)
            if (
                f in aliases or ptype not in self.ds.particle_types_raw
            ) and units not in skip_output_units:
                u = Unit(units, registry=self.ds.unit_registry)
                output_units = str(u.get_cgs_equivalent())
            else:
                output_units = units
            if (ptype, f) not in self.field_list:
                continue
            self.add_output_field(
                (ptype, f),
                sampling_type="particle",
                units=units,
                display_name=dn,
                output_units=output_units,
                take_log=False,
            )
            for alias in aliases:
                self.alias((ptype, alias), (ptype, f), units=output_units)

        ppos_fields = [f"particle_position_{ax}" for ax in "xyz"]
        pvel_fields = [f"particle_velocity_{ax}" for ax in "xyz"]
        particle_vector_functions(ptype, ppos_fields, pvel_fields, self)

        particle_deposition_functions(ptype, "particle_position", "particle_mass", self)
        standard_particle_fields(self, ptype)
        # Now we check for any leftover particle fields
        for field in sorted(self.field_list):
            if field in self:
                continue
            if not isinstance(field, tuple):
                raise RuntimeError
            if field[0] not in self.ds.particle_types:
                continue
            self.add_output_field(
                field,
                sampling_type="particle",
                units=self.ds.field_units.get(field, ""),
            )
        self.setup_smoothed_fields(ptype, num_neighbors=num_neighbors, ftype=ftype)


def _dummy_position(field, data):
    return 0.5 * np.ones_like(data[("all", "particle_position_x")])


def _dummy_velocity(field, data):
    return np.zeros_like(data[("all", "particle_velocity_x")])


def _dummy_field(field, data):
    return 0.0 * data[("chombo", "gravitational_field_x")]


fluid_field_types = ["chombo", "gas"]
particle_field_types = ["io", "all"]


class ChomboPICFieldInfo2D(ChomboPICFieldInfo3D):
    known_other_fields: KnownFieldsT = (
        ("density", (rho_units, ["density", "Density"], None)),
        (
            "potential",
            ("code_length**2 / code_time**2", ["potential", "Potential"], None),
        ),
        ("gravitational_field_x", ("code_length / code_time**2", [], None)),
        ("gravitational_field_y", ("code_length / code_time**2", [], None)),
    )
    known_particle_fields: KnownFieldsT = (
        ("particle_mass", ("code_mass", [], None)),
        ("particle_position_x", ("code_length", [], None)),
        ("particle_position_y", ("code_length", [], None)),
        ("particle_velocity_x", ("code_length / code_time", [], None)),
        ("particle_velocity_y", ("code_length / code_time", [], None)),
    )

    def __init__(self, ds, field_list):
        super().__init__(ds, field_list)

        for ftype in fluid_field_types:
            self.add_field(
                (ftype, "gravitational_field_z"),
                sampling_type="cell",
                function=_dummy_field,
                units="code_length / code_time**2",
            )

        for ptype in particle_field_types:
            self.add_field(
                (ptype, "particle_position_z"),
                sampling_type="particle",
                function=_dummy_position,
                units="code_length",
            )

            self.add_field(
                (ptype, "particle_velocity_z"),
                sampling_type="particle",
                function=_dummy_velocity,
                units="code_length / code_time",
            )


class ChomboPICFieldInfo1D(ChomboPICFieldInfo3D):
    known_other_fields: KnownFieldsT = (
        ("density", (rho_units, ["density", "Density"], None)),
        (
            "potential",
            ("code_length**2 / code_time**2", ["potential", "Potential"], None),
        ),
        ("gravitational_field_x", ("code_length / code_time**2", [], None)),
    )
    known_particle_fields: KnownFieldsT = (
        ("particle_mass", ("code_mass", [], None)),
        ("particle_position_x", ("code_length", [], None)),
        ("particle_velocity_x", ("code_length / code_time", [], None)),
    )

    def __init__(self, ds, field_list):
        super().__init__(ds, field_list)

        for ftype in fluid_field_types:
            self.add_field(
                (ftype, "gravitational_field_y"),
                sampling_type="cell",
                function=_dummy_field,
                units="code_length / code_time**2",
            )

            self.add_field(
                (ftype, "gravitational_field_z"),
                sampling_type="cell",
                function=_dummy_field,
                units="code_length / code_time**2",
            )

        for ptype in particle_field_types:
            self.add_field(
                (ptype, "particle_position_y"),
                sampling_type="particle",
                function=_dummy_position,
                units="code_length",
            )
            self.add_field(
                (ptype, "particle_position_z"),
                sampling_type="particle",
                function=_dummy_position,
                units="code_length",
            )
            self.add_field(
                (ptype, "particle_velocity_y"),
                sampling_type="particle",
                function=_dummy_velocity,
                units="code_length / code_time",
            )
            self.add_field(
                (ptype, "particle_velocity_z"),
                sampling_type="particle",
                function=_dummy_velocity,
                units="code_length / code_time",
            )


class PlutoFieldInfo(ChomboFieldInfo):
    known_other_fields: KnownFieldsT = (
        ("rho", (rho_units, ["density"], None)),
        ("prs", ("code_mass / (code_length * code_time**2)", ["pressure"], None)),
        ("vx1", (vel_units, ["velocity_x"], None)),
        ("vx2", (vel_units, ["velocity_y"], None)),
        ("vx3", (vel_units, ["velocity_z"], None)),
        ("bx1", (b_units, [], None)),
        ("bx2", (b_units, [], None)),
        ("bx3", (b_units, [], None)),
    )

    known_particle_fields = ()

    def setup_fluid_fields(self):
        from yt.fields.magnetic_field import setup_magnetic_field_aliases

        setup_magnetic_field_aliases(self, "chombo", [f"bx{ax}" for ax in [1, 2, 3]])
