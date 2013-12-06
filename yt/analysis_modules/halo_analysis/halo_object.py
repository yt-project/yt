"""
Halo Object

Author: Britton Smith <brittonsmith@gmail.com>
Affiliation: MSU
Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt-project.org/
License:
  Copyright (C) 2013 Britton Smith, Matthew Turk.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""

from yt.data_objects.data_containers import \
    YTSelectionContainer3D

class Halo(YTSelectionContainer3D):
    particles = None
    _type_name = "halo_particles"
    _skip_add = True
    def __init__(self, halo_catalog):
        self.base_source = halo_catalog.data_source
        super(Halo, self).__init__(None, halo_catalog.pf, None)

    def reset_halo(self, halo_id, particle_indices):
        self.quantities = {'halo_id':halo_id}
        self.particle_indices = particle_indices
        self.halo_id = halo_id
        self.field_data.clear()
        self.field_parameters.clear()

    def __getitem__(self, key):
        fields = self._determine_fields(key)
        if fields[0] in self.field_data:
            return self.field_data[fields[0]]
        finfo = self.pf._get_field_info(*fields[0])
        if not finfo.particle_type:
            raise NotImplementedError
        return self.base_source[fields[0]][self.particle_indices]
