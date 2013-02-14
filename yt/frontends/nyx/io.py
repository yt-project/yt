"""
Nyx data-file handling functions (basically a boxlib reader)

Author: Casey W. Stark <caseywstark@gmail.com>
Affiliation: UC Berkeley
Author: Matthew Turk <matthewturk@gmail.com>
Author: J. S. Oishi <jsoishi@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/
License:
  Copyright (C) 2011 Casey W. Stark, Matthew Turk.  All Rights Reserved.

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
from yt.utilities.io_handler import BaseIOHandler

from definitions import fab_header_pattern, nyx_particle_field_names, \
                        yt_to_nyx_fields_dict

class IOHandlerNative(BaseIOHandler):
    """ File handler that can somehow read the native boxlib format. """

    _data_style = "nyx_native"

    def modify(self, field):
        return field.swapaxes(0, 2)

    def _read_particle_field(self, grid, field):
        offset = grid._particle_offset
        filen = os.path.expanduser(grid.particle_filename)
        off = grid._particle_offset
        tr = np.zeros(grid.NumberOfParticles, dtype='float64')
        read_castro_particles(filen, off,
                              nyx_particle_field_names.index(field),
                              len(nyx_particle_field_names), tr)
        return tr

    def _read_data(self, grid, field):
        """ reads packed multiFABs output by BoxLib in "NATIVE" format. """
        if field in nyx_particle_field_names:
            return self._read_particle_field(grid, field)
        filen = os.path.expanduser(grid.filename[field])
        offset1 = grid._offset[field]
        # one field has nElements * bytesPerReal bytes and is located
        # nElements * bytesPerReal * field_index from the offset location
        bytesPerReal = grid.hierarchy._bytesPerReal

        fieldname = yt_to_nyx_fields_dict.get(field, field)
        field_index = grid.field_indexes[fieldname]
        nElements = grid.ActiveDimensions.prod()
        offset2 = int(nElements*bytesPerReal*field_index)

        dtype = grid.hierarchy._dtype
        field = np.empty(nElements, dtype=grid.hierarchy._dtype)
        read_and_seek(filen, offset1, offset2, field, nElements * bytesPerReal)
        field = field.reshape(grid.ActiveDimensions, order='F')

        # @todo: we can/should also check against the max and min in the header
        # file

        return field

