'''
AHF-specific fields



'''

#-----------------------------------------------------------------------------
# Copyright (c) 2017, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.fields.field_info_container import \
    FieldInfoContainer

m_units = 'Msun/h'
p_units = 'kpccm/h'
r_units = 'kpccm/h'
v_units = 'km/s'


class AHFHalosFieldInfo(FieldInfoContainer):
    known_other_fields = ()

    known_particle_fields = (
        ('particle_identifier', ('', ['ID'], None)),
        ('particle_mass', (m_units, ['Mvir'], 'Virial Mass')),
        ('particle_position_x', (p_units, ['Xc'], None)),
        ('particle_position_y', (p_units, ['Yc'], None)),
        ('particle_position_z', (p_units, ['Zc'], None)),
        ('particle_velocity_x', (v_units, ['VXc'], None)),
        ('particle_velocity_y', (v_units, ['VYc'], None)),
        ('particle_velocity_z', (v_units, ['VZc'], None)),
        ('virial_radius', (r_units, ['Rvir'], 'Virial Radius')),
    )
