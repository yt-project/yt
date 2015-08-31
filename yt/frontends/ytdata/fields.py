"""
YTData-specific fields




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2015, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from yt.funcs import mylog
from yt.fields.field_info_container import \
    FieldInfoContainer

m_units = "g"
p_units = "cm"
v_units = "cm / s"
r_units = "cm"

class YTDataContainerFieldInfo(FieldInfoContainer):
    known_other_fields = (
    )

    known_particle_fields = (
        ("x", (p_units, ["particle_position_x"], None)),
        ("y", (p_units, ["particle_position_y"], None)),
        ("z", (p_units, ["particle_position_z"], None)),
        ("velocity_x", (v_units, ["particle_velocity_x"], None)),
        ("velocity_y", (v_units, ["particle_velocity_y"], None)),
        ("velocity_z", (v_units, ["particle_velocity_z"], None)),
    )

    # these are extra fields to be created for the "all" particle type
    extra_union_fields = (
        (p_units, "particle_position_x"),
        (p_units, "particle_position_y"),
        (p_units, "particle_position_z"),
        (v_units, "particle_velocity_x"),
        (v_units, "particle_velocity_y"),
        (v_units, "particle_velocity_z"),
    )

class YTProjectionFieldInfo(FieldInfoContainer):
    known_other_fields = (
    )

    known_particle_fields = (
        ("px", (p_units, ["particle_position_x"], None)),
        ("py", (p_units, ["particle_position_y"], None)),
    )

class YTGridFieldInfo(FieldInfoContainer):
    known_other_fields = (
    )

    known_particle_fields = (
    )
