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
import numpy as na

class IOHandlerAthena(BaseIOHandler):
    _data_style = "athena"
    _offset_string = 'data:offsets=0'
    _data_string = 'data:datatype=0'
    _read_table_offset = None

    def _field_dict(self,fhandle):
        keys = fhandle['field_types'].keys()
        val = fhandle['field_types'].keys()
        return dict(zip(keys,val))

    def _read_field_names(self,grid):
        pass

    def _read_data_set(self,grid,field):
        f = file(grid.filename, 'rb')
        dtype, offset = grid.hierarchy._field_map[field]
        grid_ncells = na.prod(grid.ActiveDimensions)
        grid_dims = grid.ActiveDimensions
        line = f.readline()
        while True:
            splitup = line.strip().split()
            if 'CELL_DATA' in splitup:
                f.readline()
                read_table_offset = f.tell()
                del line
                break
            del line; line = f.readline()


        f.seek(read_table_offset+offset)
        if dtype == 'scalar':
            data = na.fromfile(f, dtype='>f4', count=grid_ncells).reshape(grid_dims,order='F').copy()
        if dtype == 'vector':
            data = na.fromfile(f, dtype='>f4', count=3*grid_ncells)
            if '_x' in field:
                data = data[0::3].reshape(grid_dims,order='F').copy()
            elif '_y' in field:
                data = data[1::3].reshape(grid_dims,order='F').copy()
            elif '_z' in field:
                data = data[2::3].reshape(grid_dims,order='F').copy()
        f.close()
        if grid.pf.field_ordering == 1:
            return data.T
        else:
            return data

    def _read_data_slice(self, grid, field, axis, coord):
        sl = [slice(None), slice(None), slice(None)]
        sl[axis] = slice(coord, coord + 1)
        if grid.pf.field_ordering == 1:
            sl.reverse()

        f = file(grid.filename, 'rb')
        dtype, offset = grid.hierarchy._field_map[field]
        grid_ncells = na.prod(grid.ActiveDimensions)

        line = f.readline()
        while True:
            splitup = line.strip().split()
            if 'CELL_DATA' in splitup:
                f.readline()
                read_table_offset = f.tell()
                del line
                break
            del line; line = f.readline()

        f.seek(read_table_offset+offset)
        if dtype == 'scalar':
            data = na.fromfile(f, dtype='>f4', count=grid_ncells).reshape(grid.ActiveDimensions,order='F')[sl].copy()
        if dtype == 'vector':
            data = na.fromfile(f, dtype='>f4', count=3*grid_ncells)
            if '_x' in field:
                data = data[0::3].reshape(grid.ActiveDimensions,order='F')[sl].copy()
            elif '_y' in field:
                data = data[1::3].reshape(grid.ActiveDimensions,order='F')[sl].copy()
            elif '_z' in field:
                data = data[2::3].reshape(grid.ActiveDimensions,order='F')[sl].copy()

        f.close()
        if grid.pf.field_ordering == 1:
            return data.T
        else:
            return data

