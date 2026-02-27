from yt._typing import KnownFieldsT
from yt.fields.field_info_container import FieldInfoContainer
from yt.utilities.physical_constants import (
    boltzmann_constant_cgs,
    mass_hydrogen_cgs,
)

b_units = "code_magnetic"
ra_units = "code_length / code_time**2"
rho_units = "code_density"
vel_units = "code_velocity"
pressure_units = "code_pressure"
ener_units = "code_mass * code_velocity**2"
specific_ener_units = "code_velocity**2"
ang_mom_units = "code_mass * code_velocity * code_length"


class MiniRAMSESFieldInfo(FieldInfoContainer):
    known_other_fields: KnownFieldsT = (
        ("Density", (rho_units, ["density"], None)),
        ("x-velocity", (vel_units, ["velocity_x"], None)),
        ("y-velocity", (vel_units, ["velocity_y"], None)),
        ("z-velocity", (vel_units, ["velocity_z"], None)),
        ("Pressure", (pressure_units, ["pressure"], None)),
        ("Metallicity", ("", ["metallicity"], None)),
        ("x-acceleration", (ra_units, ["acceleration_x"], None)),
        ("y-acceleration", (ra_units, ["acceleration_y"], None)),
        ("z-acceleration", (ra_units, ["acceleration_z"], None)),
        ("Potential", (specific_ener_units, ["potential"], None)),
        ("B_x", (b_units, ["magnetic_field_x"], None)),
        ("B_y", (b_units, ["magnetic_field_y"], None)),
        ("B_z", (b_units, ["magnetic_field_z"], None)),
    )
    known_particle_fields: KnownFieldsT = (
        ("particle_position_x", ("code_length", [], None)),
        ("particle_position_y", ("code_length", [], None)),
        ("particle_position_z", ("code_length", [], None)),
        ("particle_velocity_x", (vel_units, [], None)),
        ("particle_velocity_y", (vel_units, [], None)),
        ("particle_velocity_z", (vel_units, [], None)),
        ("particle_mass", ("code_mass", [], None)),
        ("particle_identity", ("", ["particle_index"], None)),
        ("particle_refinement_level", ("", [], None)),
        ("particle_birth_time", ("code_time", ["age"], None)),
        ("particle_metallicity", ("", [], None)),
    )

    def setup_particle_fields(self, ptype):
        super().setup_particle_fields(ptype)

    def setup_fluid_fields(self):
        def _temperature(data):
            rv = data["gas", "pressure"] / data["gas", "density"]
            rv *= mass_hydrogen_cgs / boltzmann_constant_cgs
            return rv

        self.add_field(
            ("gas", "temperature"),
            sampling_type="cell",
            function=_temperature,
            units=self.ds.unit_system["temperature"],
        )
