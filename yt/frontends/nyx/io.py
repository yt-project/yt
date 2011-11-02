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
import numpy as na
from yt.utilities.amr_utils import read_castro_particles
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
        tr = na.zeros(grid.NumberOfParticles, dtype='float64')
        read_castro_particles(filen, off,
                              nyx_particle_field_names.index(field),
                              len(nyx_particle_field_names), tr)
        return tr

    def _read_data_set(self, grid, field):
        """ reads packed multiFABs output by BoxLib in "NATIVE" format. """
        if field in nyx_particle_field_names:
            return self._read_particle_field(grid, field)
        filen = os.path.expanduser(grid.filename[field])
        off = grid._offset[field]
        inFile = open(filen, 'rb')
        inFile.seek(off)
        header = inFile.readline()
        header.strip()

        """
        if grid._paranoid:
            mylog.warn("Castro Native reader: Paranoid read mode.")
            header_re = re.compile(fab_header_pattern)
            bytesPerReal, endian, start, stop, centerType, nComponents = \
                headerRe.search(header).groups()

            # we will build up a dtype string, starting with endian.
            # @todo: this code is ugly.
            bytesPerReal = int(bytesPerReal)
            if bytesPerReal == int(endian[0]):
                dtype = '<'
            elif bytesPerReal == int(endian[-1]):
                dtype = '>'
            else:
                raise ValueError("FAB header is neither big nor little endian. Perhaps the file is corrupt?")

            dtype += ('f%i' % bytesPerReal)  # always a floating point

            # determine size of FAB
            start = na.array(map(int, start.split(',')))
            stop = na.array(map(int, stop.split(',')))

            gridSize = stop - start + 1

            error_count = 0
            if (start != grid.start).any():
                print "Paranoia Error: Cell_H and %s do not agree on grid start." % grid.filename
                error_count += 1
            if (stop != grid.stop).any():
                print "Paranoia Error: Cell_H and %s do not agree on grid stop." % grid.filename
                error_count += 1
            if (gridSize != grid.ActiveDimensions).any():
                print "Paranoia Error: Cell_H and %s do not agree on grid dimensions." % grid.filename
                error_count += 1
            if bytesPerReal != grid.hierarchy._bytesPerReal:
                print "Paranoia Error: Cell_H and %s do not agree on bytes per real number." % grid.filename
                error_count += 1
            if (bytesPerReal == grid.hierarchy._bytesPerReal and dtype != grid.hierarchy._dtype):
                print "Paranoia Error: Cell_H and %s do not agree on endianness." % grid.filename
                error_count += 1

            if error_count > 0:
                raise RunTimeError("Paranoia unveiled %i differences between Cell_H and %s." % (error_count, grid.filename))
        else:
        """
        start = grid.start_index
        stop = grid.stop_index
        dtype = grid.hierarchy._dtype
        bytesPerReal = grid.hierarchy._bytesPerReal

        nElements = grid.ActiveDimensions.prod()

        # one field has nElements * bytesPerReal bytes and is located
        # nElements * bytesPerReal * field_index from the offset location
        if yt_to_nyx_fields_dict.has_key(field):
            fieldname = yt_to_nyx_fields_dict[field]
        else:
            fieldname = field
        field_index = grid.field_indexes[fieldname]
        inFile.seek(int(nElements*bytesPerReal*field_index),1)
        field = na.fromfile(inFile, count=nElements, dtype=dtype)
        field = field.reshape(grid.ActiveDimensions, order='F')

        # @todo: we can/should also check against the max and min in the header
        # file

        inFile.close()
        return field

    def _read_data_slice(self, grid, field, axis, coord):
        # wishful thinking?
        sl = [slice(None), slice(None), slice(None)]
        sl[axis] = slice(coord, coord + 1)
        #sl = tuple(reversed(sl))
        return self._read_data_set(grid, field)[sl]
