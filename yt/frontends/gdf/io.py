"""
The data-file handling functions

Author: Samuel W. Skillman <samskillman@gmail.com>
Affiliation: University of Colorado at Boulder
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

    def _read_data_slice(self, grid, field, axis, coord):
        sl = [slice(None), slice(None), slice(None)]
        sl[axis] = slice(coord, coord + 1)
        if grid.pf.field_ordering == 1:
            sl.reverse()
        data = self._read_data_set(grid, field)
        if grid.pf.field_ordering == 1:
            return data.T
        else:
            return data

