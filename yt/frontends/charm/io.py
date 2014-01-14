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

from yt.utilities.io_handler import \
           BaseIOHandler

class IOHandlerCharmHDF5(BaseIOHandler):
    _data_style = "charm_hdf5"
    _offset_string = 'data:offsets=0'
    _data_string = 'data:datatype=0'

    def __init__(self, pf, *args, **kwargs):
        BaseIOHandler.__init__(self, *args, **kwargs)
        self.pf = pf
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

    def _read_particles(self, grid, name):
        num_particles = grid.pf.h.num_particles
        particles = []
        for i in np.arange(num_particles):
            particle_position_x = grid.pf.h._handle['particles/position_x'][i]
            particle_position_y = grid.pf.h._handle['particles/position_y'][i]
            particle_position_z = grid.pf.h._handle['particles/position_z'][i]
            coord = [particle_position_x, particle_position_y, particle_position_z]
            if ( (grid.LeftEdge < coord).all() and
                 (coord <= grid.RightEdge).all() ):
                particles.append(grid.pf.h._handle['particles/' + name][i])
        return np.array(particles)
