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
from yt.utilities.io_handler import \
           BaseIOHandler

from definitions import \
    yt2orionFieldsDict

class IOHandlerNative(BaseIOHandler):

    _data_style = "orion_native"

    def modify(self, field):
        return field.swapaxes(0,2)

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

        def read(line, field):
            return float(line.split(' ')[index[field]])

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

    def _read_data(self,grid,field):
        """
        reads packed multiFABs output by BoxLib in "NATIVE" format.

        """

        filen = os.path.expanduser(grid.filename[field])
        off = grid._offset[field]
        inFile = open(filen,'rb')
        inFile.seek(off)
        header = inFile.readline()
        header.strip()

        if grid._paranoid:
            mylog.warn("Orion Native reader: Paranoid read mode.")
            headerRe = re.compile(orion_FAB_header_pattern)
            bytesPerReal,endian,start,stop,centerType,nComponents = headerRe.search(header).groups()

            # we will build up a dtype string, starting with endian
            # check endianness (this code is ugly. fix?)
            bytesPerReal = int(bytesPerReal)
            if bytesPerReal == int(endian[0]):
                dtype = '<'
            elif bytesPerReal == int(endian[-1]):
                dtype = '>'
            else:
                raise ValueError("FAB header is neither big nor little endian. Perhaps the file is corrupt?")

            dtype += ('f%i'% bytesPerReal) #always a floating point

            # determine size of FAB
            start = np.array(map(int,start.split(',')))
            stop = np.array(map(int,stop.split(',')))

            gridSize = stop - start + 1

            error_count = 0
            if (start != grid.start).any():
                print "Paranoia Error: Cell_H and %s do not agree on grid start." %grid.filename
                error_count += 1
            if (stop != grid.stop).any():
                print "Paranoia Error: Cell_H and %s do not agree on grid stop." %grid.filename
                error_count += 1
            if (gridSize != grid.ActiveDimensions).any():
                print "Paranoia Error: Cell_H and %s do not agree on grid dimensions." %grid.filename
                error_count += 1
            if bytesPerReal != grid.hierarchy._bytesPerReal:
                print "Paranoia Error: Cell_H and %s do not agree on bytes per real number." %grid.filename
                error_count += 1
            if (bytesPerReal == grid.hierarchy._bytesPerReal and dtype != grid.hierarchy._dtype):
                print "Paranoia Error: Cell_H and %s do not agree on endianness." %grid.filename
                error_count += 1

            if error_count > 0:
                raise RunTimeError("Paranoia unveiled %i differences between Cell_H and %s." % (error_count, grid.filename))

        else:
            start = grid.start_index
            stop = grid.stop_index
            dtype = grid.hierarchy._dtype
            bytesPerReal = grid.hierarchy._bytesPerReal

        nElements = grid.ActiveDimensions.prod()

        # one field has nElements*bytesPerReal bytes and is located
        # nElements*bytesPerReal*field_index from the offset location
        if yt2orionFieldsDict.has_key(field):
            fieldname = yt2orionFieldsDict[field]
        else:
            fieldname = field
        field_index = grid.field_indexes[fieldname]
        inFile.seek(int(nElements*bytesPerReal*field_index),1)
        field = np.fromfile(inFile,count=nElements,dtype=dtype)
        field = field.reshape(grid.ActiveDimensions, order='F')

        # we can/should also check against the max and min in the header file

        inFile.close()
        return field

