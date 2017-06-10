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
        ("total_density", (rho_units, ["total_density"], None)),
        ("total_energy", (energy_units, ["total_energy"], None)),
        ("internal_energy", (energy_units, ["internal_energy"], None)),
    )

    known_particle_fields = (
        ("x", ("code_length", [], None)),
        ("y", ("code_length", [], None)),
        ("z", ("code_length", [], None)),
        ("vx", (vel_units, [], None)),
        ("vy", (vel_units, [], None)),
        ("vz", (vel_units, [], None)),
        ("ax", (acc_units, [], None)),
        ("ay", (acc_units, [], None)),
        ("az", (acc_units, [], None)),
        ("mass", ("code_mass", [], None)),
    )
