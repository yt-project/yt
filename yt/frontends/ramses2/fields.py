import warnings

from yt._typing import KnownFieldsT
from yt.fields.field_detector import FieldDetector
from yt.fields.field_info_container import FieldInfoContainer
from yt.utilities.physical_constants import (
    boltzmann_constant_cgs,
    mass_hydrogen_cgs,
)

b_units = "code_magnetic"
accel_units = "code_velocity / code_time"
rho_units = "code_density"
vel_units = "code_velocity"
pressure_units = "code_pressure"
specific_ener_units = "code_velocity**2"
ang_mom_units = "code_mass * code_velocity * code_length"


class RAMSES2FieldInfo(FieldInfoContainer):
    # Fields keep the names of the <prefix>_header.txt files, as in
    # utils/py/yt_ramses.py of the RAMSES II repository
    known_other_fields: KnownFieldsT = (
        ("density", (rho_units, ["density"], None)),
        ("velocity_x", (vel_units, ["velocity_x"], None)),
        ("velocity_y", (vel_units, ["velocity_y"], None)),
        ("velocity_z", (vel_units, ["velocity_z"], None)),
        ("thermal_pressure", (pressure_units, ["pressure"], None)),
        ("magnetic_field_x", (b_units, [], None)),
        ("magnetic_field_y", (b_units, [], None)),
        ("magnetic_field_z", (b_units, [], None)),
        ("metal_mass_fraction", ("", ["metallicity"], None)),
        ("xHI", ("", ["H_p0_fraction"], None)),
        ("xHII", ("", ["H_p1_fraction"], None)),
        ("xHeII", ("", ["He_p1_fraction"], None)),
        ("xHeIII", ("", ["He_p2_fraction"], None)),
        ("potential", (specific_ener_units, ["potential"], None)),
        ("accel_x", (accel_units, ["acceleration_x"], None)),
        ("accel_y", (accel_units, ["acceleration_y"], None)),
        ("accel_z", (accel_units, ["acceleration_z"], None)),
    )
    known_particle_fields: KnownFieldsT = (
        ("particle_position_x", ("code_length", [], None)),
        ("particle_position_y", ("code_length", [], None)),
        ("particle_position_z", ("code_length", [], None)),
        ("particle_velocity_x", (vel_units, [], None)),
        ("particle_velocity_y", (vel_units, [], None)),
        ("particle_velocity_z", (vel_units, [], None)),
        ("particle_acceleration_x", (accel_units, [], None)),
        ("particle_acceleration_y", (accel_units, [], None)),
        ("particle_acceleration_z", (accel_units, [], None)),
        ("particle_angular_momentum_x", (ang_mom_units, [], None)),
        ("particle_angular_momentum_y", (ang_mom_units, [], None)),
        ("particle_angular_momentum_z", (ang_mom_units, [], None)),
        ("particle_mass", ("code_mass", [], None)),
        ("particle_potential", (specific_ener_units, [], None)),
        ("particle_identity", ("", ["particle_index"], None)),
        ("particle_merging_id", ("", [], None)),
        ("particle_tracking_id", ("", [], None)),
        ("particle_level", ("", [], None)),
        ("particle_birth_time", ("code_time", ["creation_time"], None)),
        ("particle_merging_time", ("code_time", [], None)),
        ("particle_metallicity", ("", [], None)),
    )

    def setup_particle_fields(self, ptype):
        super().setup_particle_fields(ptype)

        if (ptype, "particle_birth_time") not in self.field_list:
            return

        def _age(data):
            return data.ds.current_time - data[ptype, "particle_birth_time"]

        self.add_field(
            (ptype, "age"),
            sampling_type="particle",
            function=_age,
            units=self.ds.unit_system["time"],
        )

    def setup_fluid_fields(self):
        def _temperature_over_mu(data):
            rv = data["gas", "pressure"] / data["gas", "density"]
            rv *= mass_hydrogen_cgs / boltzmann_constant_cgs
            return rv

        self.add_field(
            ("gas", "temperature_over_mu"),
            sampling_type="cell",
            function=_temperature_over_mu,
            units=self.ds.unit_system["temperature"],
        )

        def _temperature(data):
            if not isinstance(data, FieldDetector):
                warnings.warn(
                    "The mean molecular weight is not stored in RAMSES II outputs. "
                    "yt will return T/µ instead of T — this is equivalent to "
                    "assuming µ=1.0. To suppress this, derive the temperature "
                    "from temperature_over_mu with some values for mu.",
                    category=RuntimeWarning,
                    stacklevel=1,
                )
            return data["gas", "temperature_over_mu"]

        self.add_field(
            ("gas", "temperature"),
            sampling_type="cell",
            function=_temperature,
            units=self.ds.unit_system["temperature"],
        )
