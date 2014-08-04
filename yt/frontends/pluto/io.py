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

class IOHandlerPlutoHDF5(BaseIOHandler):
    _dataset_type = "pluto_hdf5"
    _offset_string = 'data:offsets=0'
    _data_string = 'data:datatype=0'

    def __init__(self, ds, *args, **kwargs):
        BaseIOHandler.__init__(self, *args, **kwargs)
        self.ds = ds
        self._handle = ds._handle

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
