import numpy as np

from yt._typing import KnownFieldsT
from yt.fields.field_info_container import FieldInfoContainer
from yt.fields.magnetic_field import setup_magnetic_field_aliases
from yt.fields.particle_fields import add_union_field
from yt.frontends.enzo_e.misc import (
    get_listed_subparam,
    get_particle_mass_correction,
    nested_dict_get,
)

rho_units = "code_mass / code_length**3"
vel_units = "code_velocity"
acc_units = "code_velocity / code_time"
energy_units = "code_velocity**2"
b_units = "code_magnetic"

NODAL_FLAGS = {
    "bfieldi_x": [1, 0, 0],
    "bfieldi_y": [0, 1, 0],
    "bfieldi_z": [0, 0, 1],
}


class EnzoEFieldInfo(FieldInfoContainer):
    known_other_fields: KnownFieldsT = (
        ("velocity_x", (vel_units, ["velocity_x"], None)),
        ("velocity_y", (vel_units, ["velocity_y"], None)),
        ("velocity_z", (vel_units, ["velocity_z"], None)),
        ("acceleration_x", (acc_units, ["acceleration_x"], None)),
        ("acceleration_y", (acc_units, ["acceleration_y"], None)),
        ("acceleration_z", (acc_units, ["acceleration_z"], None)),
        ("density", (rho_units, ["density"], None)),
        ("density_total", (rho_units, ["total_density"], None)),
        ("total_energy", (energy_units, ["specific_total_energy"], None)),
        ("internal_energy", (energy_units, ["specific_thermal_energy"], None)),
        ("bfield_x", (b_units, [], None)),
        ("bfield_y", (b_units, [], None)),
        ("bfield_z", (b_units, [], None)),
        ("bfieldi_x", (b_units, [], None)),
        ("bfieldi_y", (b_units, [], None)),
        ("bfieldi_z", (b_units, [], None)),
    )

    known_particle_fields: KnownFieldsT = (
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

    def setup_energy_field(self):
        unit_system = self.ds.unit_system
        # check if we have a field for internal energy
        has_ie_field = ("enzoe", "internal_energy") in self.field_list
        # check if we need to account for magnetic energy
        vlct_params = get_listed_subparam(self.ds.parameters, "Method", "mhd_vlct", {})
        has_magnetic = "no_bfield" != vlct_params.get("mhd_choice", "no_bfield")

        # define the ("gas", "specific_thermal_energy") field
        # - this is already done for us if the simulation used the dual-energy
        #   formalism AND ("enzoe", "internal_energy") was saved to disk
        if not (self.ds.parameters["uses_dual_energy"] and has_ie_field):

            def _tot_minus_kin(field, data):
                return (
                    data["enzoe", "total_energy"]
                    - 0.5 * data["gas", "velocity_magnitude"] ** 2.0
                )

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
        fname = "particle_mass"
        if ptype in self.ds.particle_unions:
            add_union_field(self, ptype, fname, "code_mass")
            return

        pdict = self.ds.parameters.get("Particle", None)
        if pdict is None:
            return

        constants = nested_dict_get(pdict, (ptype, "constants"), default=())
        if not constants:
            return

        # constants should be a tuple consisting of multiple tuples of (name, type, value).
        # When there is only one entry, the enclosing tuple gets stripped, so we put it back.
        if not isinstance(constants[0], tuple):
            constants = (constants,)
        names = [c[0] for c in constants]

        if "mass" not in names:
            return

        val = constants[names.index("mass")][2] * self.ds.mass_unit
        if not self.ds.index.io._particle_mass_is_mass:
            val = val * get_particle_mass_correction(self.ds)

        def _pmass(field, data):
            return val * data[ptype, "particle_ones"]

        self.add_field(
            (ptype, fname),
            function=_pmass,
            units="code_mass",
            sampling_type="particle",
        )
