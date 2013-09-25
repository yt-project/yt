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

