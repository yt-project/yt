from yt._typing import KnownFieldsT
from yt.fields.field_info_container import FieldInfoContainer
from yt.fields.magnetic_field import setup_magnetic_field_aliases
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
spin_units = "code_length * code_velocity"


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
        ("particle_spin_x", (spin_units, [], None)),
        ("particle_spin_y", (spin_units, [], None)),
        ("particle_spin_z", (spin_units, [], None)),
        ("particle_mass", ("code_mass", [], None)),
        ("particle_potential", (specific_ener_units, [], None)),
        ("particle_identity", ("", ["particle_index"], None)),
        ("particle_merging_id", ("", [], None)),
        ("particle_tracking_id", ("", [], None)),
        ("particle_level", ("", [], None)),
        ("particle_birth_time", ("code_time", ["creation_time"], None)),
        ("particle_merging_time", ("code_time", [], None)),
        ("particle_metallicity", ("", [], None)),
        ("particle_internal_energy", (specific_ener_units, [], None)),
        ("particle_accretion_step", ("", [], None)),
    )

    def setup_particle_fields(self, ptype):
        super().setup_particle_fields(ptype)

        if (ptype, "particle_birth_time") in self.field_list:

            def _age(data):
                return data.ds.current_time - data[ptype, "particle_birth_time"]

            self.add_field(
                (ptype, "age"),
                sampling_type="particle",
                function=_age,
                units=self.ds.unit_system["time"],
            )

    def setup_fluid_fields(self):
        setup_magnetic_field_aliases(
            self, "ramses", [f"magnetic_field_{ax}" for ax in "xyz"]
        )

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
