"""
Fields specific to Enzo-P



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.fields.field_info_container import \
    FieldInfoContainer
from yt.fields.particle_fields import \
    add_union_field
from yt.frontends.enzo_p.misc import \
    nested_dict_get

rho_units = "code_mass / code_length**3"
vel_units = "code_velocity"
acc_units = "code_velocity / code_time"
energy_units = "code_velocity**2"

known_species_names = {
}

class EnzoPFieldInfo(FieldInfoContainer):
    known_other_fields = (
        ("velocity_x", (vel_units, ["velocity_x"], None)),
        ("velocity_y", (vel_units, ["velocity_y"], None)),
        ("velocity_z", (vel_units, ["velocity_z"], None)),
        ("acceleration_x", (acc_units, ["acceleration_x"], None)),
        ("acceleration_y", (acc_units, ["acceleration_y"], None)),
        ("acceleration_z", (acc_units, ["acceleration_z"], None)),
        ("density", (rho_units, ["density"], None)),
        ("density_total", (rho_units, ["total_density"], None)),
        ("total_energy", (energy_units, ["total_energy"], None)),
        ("internal_energy", (energy_units, ["internal_energy"], None)),
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

    def __init__(self, ds, field_list, slice_info = None):
        super(EnzoPFieldInfo, self).__init__(
            ds, field_list, slice_info=slice_info)

    def setup_particle_fields(self, ptype, ftype='gas', num_neighbors=64):
        super(EnzoPFieldInfo, self).setup_particle_fields(
            ptype, ftype=ftype, num_neighbors=num_neighbors)
        self.setup_particle_mass_field(ptype)

    def setup_particle_mass_field(self, ptype):
        name = "particle_mass"
        if ptype in self.ds.particle_unions:
            add_union_field(self, ptype, name, "code_mass")
            return

        constants = nested_dict_get(
            self.ds.parameters, ("Particle", ptype, "constants"),
            default=[])
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
                val /= self.ds.domain_dimensions.prod()

            def _pmass(field, data):
                return val * data[ptype, "particle_ones"]
            self.add_field((ptype, name),
                            function=_pmass, units="code_mass",
                            sampling_type="particle")
