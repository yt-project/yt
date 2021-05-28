import numpy as np

from yt.fields.field_info_container import FieldInfoContainer
from yt.fields.magnetic_field import setup_magnetic_field_aliases
from yt.fields.particle_fields import add_union_field
from yt.frontends.enzo_e.misc import nested_dict_get

rho_units = "code_mass / code_length**3"
vel_units = "code_velocity"
acc_units = "code_velocity / code_time"
energy_units = "code_velocity**2"
b_units = "code_magnetic"

known_species_names = {}

NODAL_FLAGS = {
    "bfieldi_x": [1, 0, 0],
    "bfieldi_y": [0, 1, 0],
    "bfieldi_z": [0, 0, 1],
}


class EnzoEFieldInfo(FieldInfoContainer):
    known_other_fields = (
        ("velocity_x", (vel_units, ["velocity_x"], None)),
        ("velocity_y", (vel_units, ["velocity_y"], None)),
        ("velocity_z", (vel_units, ["velocity_z"], None)),
        ("acceleration_x", (acc_units, ["acceleration_x"], None)),
        ("acceleration_y", (acc_units, ["acceleration_y"], None)),
        ("acceleration_z", (acc_units, ["acceleration_z"], None)),
        ("density", (rho_units, ["density"], None)),
        ("density_total", (rho_units, ["total_density"], None)),
        ("total_energy", (energy_units, ["specific_total_energy"], None)),
        ("internal_energy", (energy_units, ["specific_internal_energy"], None)),
        ("bfield_x", (b_units, [], None)),
        ("bfield_y", (b_units, [], None)),
        ("bfield_z", (b_units, [], None)),
        ("bfieldi_x", (b_units, [], None)),
        ("bfieldi_y", (b_units, [], None)),
        ("bfieldi_z", (b_units, [], None)),
    )

    known_particle_fields = (
        ("x", ("code_length", ["particle_position_x"], None)),
        ("y", ("code_length", ["particle_position_y"], None)),
        ("z", ("code_length", ["particle_position_z"], None)),
        ("vx", (vel_units, ["particle_velocity_x"], None)),
        ("vy", (vel_units, ["particle_velocity_y"], None)),
        ("vz", (vel_units, ["particle_velocity_z"], None)),
        ("ax", (acc_units, ["particle_acceleration_x"], None)),
        ("ay", (acc_units, ["particle_acceleration_y"], None)),
        ("az", (acc_units, ["particle_acceleration_z"], None)),
        ("mass", ("code_mass", ["particle_mass"], None)),
    )

    def __init__(self, ds, field_list, slice_info=None):
        super().__init__(ds, field_list, slice_info=slice_info)

        # setup nodal flag information
        for field, arr in NODAL_FLAGS.items():
            if ("enzoe", field) in self:
                finfo = self["enzoe", field]
                finfo.nodal_flag = np.array(arr)

    def setup_fluid_fields(self):
        self.setup_energy_field()
        setup_magnetic_field_aliases(self, "enzoe", [f"bfield_{ax}" for ax in "xyz"])
        self.alias(
            ("gas", "total_energy"),
            ("gas", "specific_total_energy"),
            deprecate=("4.0.0", "4.1.0"),
        )

    def uses_dual_energy_formalism(self):
        for name in ["ppm", "mhd_vlct"]:
            param = nested_dict_get(self.ds.parameters, ("Method", name), None)
            if (param is not None) and param.get("dual_energy", False):
                return True
        return False

    def setup_energy_field(self):
        unit_system = self.ds.unit_system
        # check if we need to include magnetic energy
        has_magnetic = "bfield_x" in self.ds.parameters["Field"]["list"]

        # identify if the dual energy formalism is in use:
        if self.uses_dual_energy_formalism():
            self.alias(
                ("gas", "specific_thermal_energy"),
                ("enzop", "specific_internal_energy"),
            )
        else:

            def _tot_minus_kin(field, data):
                return (
                    data["enzoe", "total_energy"]
                    - 0.5 * data["gas", "velocity_magnitude"] ** 2.0
                )

            # check if we need to include magnetic energy
            has_magnetic = "bfield_x" in self.ds.parameters["Field"]["list"]

            if not has_magnetic:
                # thermal energy = total energy - kinetic energy
                self.add_field(
                    ("gas", "specific_thermal_energy"),
                    sampling_type="cell",
                    function=_tot_minus_kin,
                    units=unit_system["specific_energy"],
                )
            else:
                # thermal energy = total energy - kinetic energy - magnetic energy
                def _sub_b(field, data):
                    return (
                        _tot_minus_kin(field, data)
                        - data["gas", "magnetic_energy_density"]
                        / data["gas", "density"]
                    )

                self.add_field(
                    ("gas", "specific_thermal_energy"),
                    sampling_type="cell",
                    function=_sub_b,
                    units=unit_system["specific_energy"],
                )

    def setup_particle_fields(self, ptype, ftype="gas", num_neighbors=64):
        super().setup_particle_fields(ptype, ftype=ftype, num_neighbors=num_neighbors)
        self.setup_particle_mass_field(ptype)

    def setup_particle_mass_field(self, ptype):
        name = "particle_mass"
        if ptype in self.ds.particle_unions:
            add_union_field(self, ptype, name, "code_mass")
            return

        constants = nested_dict_get(
            self.ds.parameters, ("Particle", ptype, "constants"), default=[]
        )
        if not constants:
            names = []
        else:
            if not isinstance(constants[0], tuple):
                constants = (constants,)
            names = [c[0] for c in constants]

        if "mass" in names:
            val = constants[names.index("mass")][2]
            val = self.ds.quan(val, self.ds.mass_unit)
            if self.ds.cosmological_simulation:
                val = val / self.ds.domain_dimensions.prod()

            def _pmass(field, data):
                return val * data[ptype, "particle_ones"]

            self.add_field(
                (ptype, name),
                function=_pmass,
                units="code_mass",
                sampling_type="particle",
            )
