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

        particle_field_index = {'position_x': 0,
                                'position_y': 1,
                                'position_z': 2,
                                'velocity_x': 3,
                                'velocity_y': 4,
                                'velocity_z': 5,
                                'acceleration_x': 6,
                                'acceleration_y': 7,
                                'acceleration_z': 8,
                                'mass': 9}

        field_index = particle_field_index[name]

        particles_per_cell = self._handle['level_0/particles:offsets'].value
        items_per_particle = len(particle_field_index)

        # compute global offset position
        offsets = items_per_particle * np.cumsum(particles_per_cell)
        offsets = np.append(np.array([0]), offsets)
        offsets = np.array(offsets, dtype=np.int64)
        
    
        data = self._handle['level_0/particles:data'][offsets[grid.id]:offsets[grid.id + 1]]
        return data[field_index::items_per_particle]

class IOHandlerCharm2DHDF5(IOHandlerCharmHDF5):
    _data_style = "charm2d_hdf5"
    _offset_string = 'data:offsets=0'
    _data_string = 'data:datatype=0'

    def __init__(self, pf, *args, **kwargs):
        BaseIOHandler.__init__(self, *args, **kwargs)
        self.pf = pf
        self._handle = pf._handle

    def _read_particles(self, grid, name):

        particle_field_index = {'position_x': 0,
                                'position_y': 1,
                                'velocity_x': 2,
                                'velocity_y': 3,
                                'acceleration_x': 4,
                                'acceleration_y': 5,
                                'mass': 6}

        field_index = particle_field_index[name]

        particles_per_cell = self._handle['level_0/particles:offsets'].value
        items_per_particle = len(particle_field_index)

        # compute global offset position
        offsets = items_per_particle * np.cumsum(particles_per_cell)
        offsets = np.append(np.array([0]), offsets)
        offsets = np.array(offsets, dtype=np.int64)
           
        data = self._handle['level_0/particles:data'][offsets[grid.id]:offsets[grid.id + 1]]
        return data[field_index::items_per_particle]

class IOHandlerCharm1DHDF5(IOHandlerCharmHDF5):
    _data_style = "charm1d_hdf5"
    _offset_string = 'data:offsets=0'
    _data_string = 'data:datatype=0'

    def __init__(self, pf, *args, **kwargs):
        BaseIOHandler.__init__(self, *args, **kwargs)
        self.pf = pf
        self._handle = pf._handle

    def _read_particles(self, grid, name):

        particle_field_index = {'position_x': 0,
                                'velocity_x': 1,
                                'acceleration_x': 2,
                                'mass': 3}

        field_index = particle_field_index[name]

        particles_per_cell = self._handle['level_0/particles:offsets'].value
        items_per_particle = len(particle_field_index)

        # compute global offset position
        offsets = items_per_particle * np.cumsum(particles_per_cell)
        offsets = np.append(np.array([0]), offsets)
        offsets = np.array(offsets, dtype=np.int64)
           
        data = self._handle['level_0/particles:data'][offsets[grid.id]:offsets[grid.id + 1]]
        return data[field_index::items_per_particle]