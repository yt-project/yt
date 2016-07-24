"""
SDF-specific fields




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


class SDFFieldInfo(FieldInfoContainer):
    known_other_fields = ()

    known_particle_fields = ()
    _mass_field = None

    def __init__(self, ds, field_list):

        if 'mass' in field_list:
            self.known_particle_fields.append(("mass", "code_mass",
                                               ["particle_mass"], None))
        possible_masses = ['mass', 'm200b', 'mvir']
        mnf = 'mass'
        for mn in possible_masses:
            if mn in ds.sdf_container.keys():
                mnf = self._mass_field = mn
                break

        idf = ds._field_map.get("particle_index", 'ident')
        xf = ds._field_map.get("particle_position_x", 'x')
        yf = ds._field_map.get("particle_position_y", 'y')
        zf = ds._field_map.get("particle_position_z", 'z')
        vxf = ds._field_map.get("particle_velocity_x", 'vx')
        vyf = ds._field_map.get("particle_velocity_z", 'vy')
        vzf = ds._field_map.get("particle_velocity_z", 'vz')

        self.known_particle_fields = (
            (idf, ('dimensionless', ['particle_index'], None)),
            (xf,  ('code_length', ['particle_position_x'], None)),
            (yf,  ('code_length', ['particle_position_y'], None)),
            (zf,  ('code_length', ['particle_position_z'], None)),
            (vxf, ('code_velocity', ['particle_velocity_x'], None)),
            (vyf, ('code_velocity', ['particle_velocity_y'], None)),
            (vzf, ('code_velocity', ['particle_velocity_z'], None)),
            (mnf, ('code_mass', ['particle_mass'], None)),
        )
        super(SDFFieldInfo, self).__init__(ds, field_list)


