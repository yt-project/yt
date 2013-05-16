"""
Orion data-file handling functions

Author: Matthew Turk <matthewturk@gmail.com>
Author: J. S. Oishi <jsoishi@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/
License:
  Copyright (C) 2007-2011 Matthew Turk.  All Rights Reserved.

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

import os
import numpy as np
from yt.utilities.lib import read_castro_particles, read_and_seek
from yt.utilities.io_handler import \
           BaseIOHandler

class IOHandlerBoxlib(BaseIOHandler):

    _data_style = "boxlib_native"

    def modify(self, field):
        return field.swapaxes(0,2)

        def read(line, field):
            return float(line.split(' ')[index[field]])

    def _read_data(self, grid, field):
        """
        reads packed multiFABs output by BoxLib in "NATIVE" format.

        """

        filen = os.path.expanduser(grid.filename)
        offset1 = grid._offset
        # one field has nElements * bytesPerReal bytes and is located
        # nElements * bytesPerReal * field_index from the offset location
        bytesPerReal = grid.hierarchy._bytesPerReal

        field_index = grid.hierarchy.field_indexes[field]
        nElements = grid.ActiveDimensions.prod()
        offset2 = int(nElements*bytesPerReal*field_index)

        dtype = grid.hierarchy._dtype
        field = np.empty(nElements, dtype=grid.hierarchy._dtype)
        read_and_seek(filen, offset1, offset2, field, nElements * bytesPerReal)
        field = field.reshape(grid.ActiveDimensions, order='F')

        return field

class IOHandlerOrion(IOHandlerBoxlib):
    _data_style = "orion_native"

    def _read_particles(self, grid, field): 
        """
        parses the Orion Star Particle text files
        
        """
        index = {'particle_mass': 0,
                 'particle_position_x': 1,
                 'particle_position_y': 2,
                 'particle_position_z': 3,
                 'particle_momentum_x': 4,
                 'particle_momentum_y': 5,
                 'particle_momentum_z': 6,
                 'particle_angmomen_x': 7,
                 'particle_angmomen_y': 8,
                 'particle_angmomen_z': 9,
                 'particle_mlast': 10,
                 'particle_mdeut': 11,
                 'particle_n': 12,
                 'particle_mdot': 13,
                 'particle_burnstate': 14,
                 'particle_id': 15}

        fn = grid.pf.fullplotdir + "/StarParticles"
        with open(fn, 'r') as f:
            lines = f.readlines()
            particles = []
            for line in lines[1:]:
                if grid.NumberOfParticles > 0:
                    coord = read(line, "particle_position_x"), \
                            read(line, "particle_position_y"), \
                            read(line, "particle_position_z") 
                    if ( (grid.LeftEdge < coord).all() and 
                         (coord <= grid.RightEdge).all() ):
                        particles.append(read(line, field))
        return np.array(particles)

class IOHandlerCastro(IOHandlerBoxlib):
    _data_style = "castro_native"

    def _read_particle_field(self, grid, field):
        offset = grid._particle_offset
        filen = os.path.expanduser(grid.particle_filename)
        off = grid._particle_offset
        tr = np.zeros(grid.NumberOfParticles, dtype='float64')
        read_castro_particles(filen, off,
            castro_particle_field_names.index(field),
            len(castro_particle_field_names),
            tr)
        return tr

