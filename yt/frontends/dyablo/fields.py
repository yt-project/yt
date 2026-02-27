"""
Field definitions for Dyablo frontend.
"""

import unyt

from yt._typing import KnownFieldsT
from yt.fields.field_info_container import FieldInfoContainer

# Unit definitions
rho_units = "code_density"
mom_units = "code_density * code_velocity"
en_units = "code_density * code_velocity**2"
vel_units = "code_velocity"


class DyabloFieldInfo(FieldInfoContainer):
    """Field info container for Dyablo code."""

    known_other_fields: KnownFieldsT = (
        ("rho", (rho_units, ["density"], None)),
        ("rho_vx", (mom_units, ["density_times_velocity_x"], None)),
        ("rho_vy", (mom_units, ["density_times_velocity_y"], None)),
        ("rho_vz", (mom_units, ["density_times_velocity_z"], None)),
        ("e_tot", (en_units, ["total_energy_density"], None)),
        ("metallicity", ("dimensionless", ["metallicity"], None)),
    )

    known_particle_fields: KnownFieldsT = (
        ("particle_position_x", ("code_length", [], None)),
        ("particle_position_y", ("code_length", [], None)),
        ("particle_position_z", ("code_length", [], None)),
        ("particle_identity", ("", ["particle_index"], None)),
        ("particle_vx", (vel_units, ["particle_velocity_x"], None)),
        ("particle_vy", (vel_units, ["particle_velocity_y"], None)),
        ("particle_vz", (vel_units, ["particle_velocity_z"], None)),
        ("particle_mass", ("code_mass", [], None)),
        ("particle_birth_time", ("code_time", [], None)),
    )

    def setup_fluid_fields(self):
        """Set up derived fluid fields from conserved quantities."""

        # Create velocity components from momentum
        def _velocity_x(data):
            return data["dyablo", "rho_vx"] / data["dyablo", "rho"]

        def _velocity_y(data):
            return data["dyablo", "rho_vy"] / data["dyablo", "rho"]

        def _velocity_z(data):
            return data["dyablo", "rho_vz"] / data["dyablo", "rho"]

        self.add_field(
            ("gas", "velocity_x"),
            sampling_type="cell",
            function=_velocity_x,
            units=self.ds.unit_system["velocity"],
        )
        self.add_field(
            ("gas", "velocity_y"),
            sampling_type="cell",
            function=_velocity_y,
            units=self.ds.unit_system["velocity"],
        )
        self.add_field(
            ("gas", "velocity_z"),
            sampling_type="cell",
            function=_velocity_z,
            units=self.ds.unit_system["velocity"],
        )

        # Create velocity magnitude
        def _velocity_magnitude(data):
            return (
                data["gas", "velocity_x"] ** 2
                + data["gas", "velocity_y"] ** 2
                + data["gas", "velocity_z"] ** 2
            ) ** 0.5

        self.add_field(
            ("gas", "velocity_magnitude"),
            sampling_type="cell",
            function=_velocity_magnitude,
            units=self.ds.unit_system["velocity"],
        )

        # Create pressure
        gamma = self.ds.parameters.get("gamma", 5.0 / 3.0)

        def _pressure(data):
            kinetic_energy = 0.5 * (
                data["gas", "velocity_x"] ** 2
                + data["gas", "velocity_y"] ** 2
                + data["gas", "velocity_z"] ** 2
            )
            internal_energy = (
                data["dyablo", "e_tot"] / data["dyablo", "rho"] - kinetic_energy
            )
            return (gamma - 1.0) * data["dyablo", "rho"] * internal_energy

        self.add_field(
            ("gas", "pressure"),
            sampling_type="cell",
            function=_pressure,
            units=self.ds.unit_system["pressure"],
        )

        # Create temperature
        def _temperature_over_mu(data):
            rv = data["gas", "pressure"] / data["gas", "density"]
            return rv * unyt.mp / unyt.kb

        self.add_field(
            ("gas", "temperature_over_mu"),
            sampling_type="cell",
            function=_temperature_over_mu,
            units=self.ds.unit_system["temperature"],
        )

    def setup_particle_fields(self, ptype):
        """Set up particle fields."""
        super().setup_particle_fields(ptype)
