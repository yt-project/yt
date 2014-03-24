"""
The data-file handling functions



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------
import h5py
import os
import re
import numpy as np
from yt.utilities.logger import ytLogger as mylog

from yt.utilities.io_handler import \
           BaseIOHandler

class IOHandlerChomboHDF5(BaseIOHandler):
    _dataset_type = "chombo_hdf5"
    _offset_string = 'data:offsets=0'
    _data_string = 'data:datatype=0'

    def __init__(self, pf, *args, **kwargs):
        BaseIOHandler.__init__(self, pf)
        self._handle = pf._handle

    _field_dict = None
    @property
    def field_dict(self):
        if self._field_dict is not None:
            return self._field_dict
        ncomp = int(self._handle['/'].attrs['num_components'])
        temp =  self._handle['/'].attrs.items()[-ncomp:]
        val, keys = zip(*temp)
        val = [int(re.match('component_(\d+)',v).groups()[0]) for v in val]
        self._field_dict = dict(zip(keys,val))
        return self._field_dict
        
    def _read_field_names(self,grid):
        ncomp = int(self._handle['/'].attrs['num_components'])

        fns = [c[1] for c in f['/'].attrs.items()[-ncomp-1:-1]]
    
    def _read_data(self,grid,field):

        lstring = 'level_%i' % grid.Level
        lev = self._handle[lstring]
        dims = grid.ActiveDimensions
        boxsize = dims.prod()
        
        grid_offset = lev[self._offset_string][grid._level_id]
        start = grid_offset+self.field_dict[field]*boxsize
        stop = start + boxsize
        data = lev[self._data_string][start:stop]
        
        return data.reshape(dims, order='F')

    def _read_fluid_selection(self, chunks, selector, fields, size):
        rv = {}
        chunks = list(chunks)
        fields.sort(key=lambda a: self.field_dict[a[1]])
        if selector.__class__.__name__ == "GridSelector":
            if not (len(chunks) == len(chunks[0].objs) == 1):
                raise RuntimeError
            grid = chunks[0].objs[0]
            lstring = 'level_%i' % grid.Level
            lev = self._handle[lstring]
            grid_offset = lev[self._offset_string][grid._level_id]
            boxsize = grid.ActiveDimensions.prod()
            for ftype, fname in fields:
                start = grid_offset+self.field_dict[fname]*boxsize
                stop = start + boxsize
                data = lev[self._data_string][start:stop]
                rv[ftype, fname] = data.reshape(grid.ActiveDimensions,
                                        order='F')
            return rv
        if size is None:
            size = sum((g.count(selector) for chunk in chunks
                        for g in chunk.objs))
        for field in fields:
            ftype, fname = field
            fsize = size
            rv[field] = np.empty(fsize, dtype="float64")
        ng = sum(len(c.objs) for c in chunks)
        mylog.debug("Reading %s cells of %s fields in %s grids",
                   size, [f2 for f1, f2 in fields], ng)
        ind = 0
        for chunk in chunks:
            for g in chunk.objs:
                lstring = 'level_%i' % g.Level
                lev = self._handle[lstring]
                grid_offset = lev[self._offset_string][g._level_id]
                boxsize = g.ActiveDimensions.prod()
                nd = 0
                for field in fields:
                    start = grid_offset+self.field_dict[fname]*boxsize
                    stop = start + boxsize
                    data = lev[self._data_string][start:stop]
                    data = data.reshape(g.ActiveDimensions, order='F')
                    nd = g.select(selector, data, rv[field], ind) # caches
                ind += nd
        return rv

    def _read_particles(self, grid, field):
        """
        parses the Orion Star Particle text files
             
        """

        fn = grid.pf.fullplotdir[:-4] + "sink"

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

        fn = grid.pf.fullplotdir[:-4] + "sink"
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
