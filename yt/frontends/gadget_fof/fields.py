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
p_units = "code_length"
v_units = "code_velocity"

class GadgetFOFFieldInfo(FieldInfoContainer):
    known_other_fields = (
    )

    known_particle_fields = (
        ("GroupPos_0", (p_units, ["particle_position_x"], None)),
        ("GroupPos_1", (p_units, ["particle_position_y"], None)),
        ("GroupPos_2", (p_units, ["particle_position_z"], None)),
        ("GroupVel_0", (v_units, ["particle_velocity_x"], None)),
        ("GroupVel_1", (v_units, ["particle_velocity_y"], None)),
        ("GroupVel_2", (v_units, ["particle_velocity_z"], None)),
        ("GroupMass",  (m_units, ["particle_mass"], None)),
        ("GroupLen",   ("",      ["particle_number"], None)),
)
