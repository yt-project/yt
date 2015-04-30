"""
GadgetFOF-specific fields




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2015, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.funcs import mylog
from yt.fields.field_info_container import \
    FieldInfoContainer
from yt.units.yt_array import \
    YTArray

m_units = "code_mass"
mdot_units = "code_mass / code_time"
p_units = "Mpccm/h"
v_units = "1e5 * cmcm / s"

class GadgetFOFFieldInfo(FieldInfoContainer):
    known_other_fields = (
    )

    known_particle_fields = (
        ("Pos_0", (p_units, ["particle_position_x"], None)),
        ("Pos_1", (p_units, ["particle_position_y"], None)),
        ("Pos_2", (p_units, ["particle_position_z"], None)),
        ("Vel_0", (v_units, ["particle_velocity_x"], None)),
        ("Vel_1", (v_units, ["particle_velocity_y"], None)),
        ("Vel_2", (v_units, ["particle_velocity_z"], None)),
        ("Mass", (m_units, ["particle_mass"], None)),
)
