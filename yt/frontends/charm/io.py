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
        self._particle_field_index = {'position_x': 0,
                                      'position_y': 1,
                                      'position_z': 2,
                                      'velocity_x': 3,
                                      'velocity_y': 4,
                                      'velocity_z': 5,
                                      'acceleration_x': 6,
                                      'acceleration_y': 7,
                                      'acceleration_z': 8,
                                      'mass': 9}

    _field_dict = None
    @property
    def field_dict(self):
        if self._field_dict is not None:
            return self._field_dict
        field_dict = {}
        for key, val in self._handle.attrs.items():
            if key.startswith('component_'):
                comp_number = int(re.match('component_(\d)', key).groups()[0])
                field_dict[val] = comp_number
        self._field_dict = field_dict
        return self._field_dict
        
    def _read_field_names(self, grid):
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

        field_index = self._particle_field_index[name]
        lev = 'level_%s' % grid.Level

        particles_per_grid = self._handle[lev]['particles:offsets'].value
        items_per_particle = len(self._particle_field_index)

        # compute global offset position
        offsets = items_per_particle * np.cumsum(particles_per_grid)
        offsets = np.append(np.array([0]), offsets)
        offsets = np.array(offsets, dtype=np.int64)

        # convert between the global grid id and the id on this level            
        grid_levels = np.array([g.Level for g in self.pf.h.grids])
        grid_ids    = np.array([g.id    for g in self.pf.h.grids])
        grid_level_offset = grid_ids[np.where(grid_levels == grid.Level)[0][0]]
        lo = grid.id - grid_level_offset
        hi = lo + 1

        data = self._handle[lev]['particles:data'][offsets[lo]:offsets[hi]]
        return data[field_index::items_per_particle]

class IOHandlerCharm2DHDF5(IOHandlerCharmHDF5):
    _data_style = "charm2d_hdf5"
    _offset_string = 'data:offsets=0'
    _data_string = 'data:datatype=0'

    def __init__(self, pf, *args, **kwargs):
        BaseIOHandler.__init__(self, *args, **kwargs)
        self.pf = pf
        self._handle = pf._handle
        self._particle_field_index = {'position_x': 0,
                                      'position_y': 1,
                                      'velocity_x': 2,
                                      'velocity_y': 3,
                                      'acceleration_x': 4,
                                      'acceleration_y': 5,
                                      'mass': 6}


class IOHandlerCharm1DHDF5(IOHandlerCharmHDF5):
    _data_style = "charm1d_hdf5"
    _offset_string = 'data:offsets=0'
    _data_string = 'data:datatype=0'

    def __init__(self, pf, *args, **kwargs):
        BaseIOHandler.__init__(self, *args, **kwargs)
        self.pf = pf
        self._handle = pf._handle
        self._particle_field_index = {'position_x': 0,
                                      'velocity_x': 1,
                                      'acceleration_x': 2,
                                      'mass': 3}