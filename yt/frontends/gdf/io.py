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
from yt.utilities.io_handler import \
           BaseIOHandler
import h5py

class IOHandlerGDFHDF5(BaseIOHandler):
    _data_style = "grid_data_format"
    _offset_string = 'data:offsets=0'
    _data_string = 'data:datatype=0'

    def _field_dict(self,fhandle):
        keys = fhandle['field_types'].keys()
        val = fhandle['field_types'].keys()
        return dict(zip(keys,val))

    def _read_field_names(self,grid):
        fhandle = h5py.File(grid.filename,'r')
        names = fhandle['field_types'].keys()
        fhandle.close()
        return names

    def _read_data(self,grid,field):
        fhandle = h5py.File(grid.hierarchy.hierarchy_filename,'r')
        data = (fhandle['/data/grid_%010i/'%grid.id+field][:]).copy()
        fhandle.close()
        if grid.pf.field_ordering == 1:
            return data.T
        else:
            return data
