"""
Orion data-file handling functions



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import os
import numpy as np
from yt.utilities.lib import read_castro_particles, read_and_seek
from yt.utilities.io_handler import \
           BaseIOHandler

class IOHandlerBoxlib(BaseIOHandler):

    _data_style = "boxlib_native"

    def modify(self, field):
        return field.swapaxes(0,2)

    def _read_data(self, grid, field):
        """
        reads packed multiFABs output by BoxLib in "NATIVE" format.

        """

        filen = os.path.expanduser(grid.filename)
        offset1 = grid._offset
        dtype = grid.hierarchy._dtype
        bpr = grid.hierarchy._dtype.itemsize
        ne = grid.ActiveDimensions.prod()
        field_index = grid.hierarchy.field_indexes[field]
        offset2 = int(ne*bpr*field_index)

        field = np.empty(grid.ActiveDimensions, dtype=dtype, order='F')
        read_and_seek(filen, offset1, offset2, field, ne * bpr)

        return field

class IOHandlerOrion(IOHandlerBoxlib):
    _data_style = "orion_native"

    def _read_particles(self, grid, field): 
        """
        parses the Orion Star Particle text files
        
        """

        fn = grid.pf.fullplotdir + "/StarParticles"
        if not os.path.exists(fn):
            fn = grid.pf.fullplotdir + "/SinkParticles"

        # Figure out the format of the particle file
        with open(fn, 'r') as f:
            lines = f.readlines()
        line = lines[1]
        
        # The basic fields that all sink particles have
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
                 'particle_id': -1}

        if len(line.strip().split()) == 11:
            # these are vanilla sinks, do nothing
            pass  

        elif len(line.strip().split()) == 17:
            # these are old-style stars, add stellar model parameters
            index['particle_mlast']     = 10
            index['particle_r']         = 11
            index['particle_mdeut']     = 12
            index['particle_n']         = 13
            index['particle_mdot']      = 14,
            index['particle_burnstate'] = 15

        elif len(line.strip().split()) == 18:
            # these are the newer style, add luminosity as well
            index['particle_mlast']     = 10
            index['particle_r']         = 11
            index['particle_mdeut']     = 12
            index['particle_n']         = 13
            index['particle_mdot']      = 14,
            index['particle_burnstate'] = 15,
            index['particle_luminosity']= 16

        else:
            # give a warning if none of the above apply:
            mylog.warning('Warning - could not figure out particle output file')
            mylog.warning('These results could be nonsense!')

        def read(line, field):
            return float(line.strip().split(' ')[index[field]])

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

