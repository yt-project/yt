"""
ytree-specific fields




"""

#-----------------------------------------------------------------------------
# Copyright (c) yt Development Team. All rights reserved.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.fields.field_info_container import \
    FieldInfoContainer

p_units = 'unitary'
v_units = 'km/s'

class YTreeFieldInfo(FieldInfoContainer):
    known_other_fields = (
    )

    known_particle_fields = (
        ("position_x", (p_units, ('halos', 'particle_position_x'), None)),
        ("position_y", (p_units, ('halos', 'particle_position_y'), None)),
        ("position_z", (p_units, ('halos', 'particle_position_z'), None)),
        ("velocity_x", (v_units, ('halos', 'particle_velocity_x'), None)),
        ("velocity_y", (v_units, ('halos', 'particle_velocity_y'), None)),
        ("velocity_z", (v_units, ('halos', 'particle_velocity_z'), None)),
        ("mass", ('Msun', ('halos', 'particle_mass'), None)),
    )
