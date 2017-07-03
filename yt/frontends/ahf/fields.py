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
        ('ID', ('', ['particle_identifier'], None)),
        ('Mvir', (m_units, ['particle_mass'], 'Virial Mass')),
        ('Xc', (p_units, ['particle_position_x'], None)),
        ('Yc', (p_units, ['particle_position_y'], None)),
        ('Zc', (p_units, ['particle_position_z'], None)),
        ('VXc', (v_units, ['particle_velocity_x'], None)),
        ('VYc', (v_units, ['particle_velocity_y'], None)),
        ('VZc', (v_units, ['particle_velocity_z'], None)),
        ('Rvir', (r_units, ['virial_radius'], 'Virial Radius')),
    )
