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
import numpy as np

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

    def _read_data(self,grid,field):
        f = file(grid.filename, 'rb')
        dtype, offsetr = grid.hierarchy._field_map[field]
        grid_ncells = np.prod(grid.ActiveDimensions)
        grid_dims = grid.ActiveDimensions
        grid0_ncells = np.prod(grid.hierarchy.grid_dimensions[0,:])
        read_table_offset = get_read_table_offset(f)
        if grid_ncells != grid0_ncells:
            offset = offsetr + ((grid_ncells-grid0_ncells) * (offsetr//grid0_ncells))
        if grid_ncells == grid0_ncells:
            offset = offsetr
        f.seek(read_table_offset+offset)
        if dtype == 'scalar':
            data = np.fromfile(f, dtype='>f4',
                    count=grid_ncells).reshape(grid_dims,order='F').copy()
        if dtype == 'vector':
            data = np.fromfile(f, dtype='>f4', count=3*grid_ncells)
            if '_x' in field:
                data = data[0::3].reshape(grid_dims,order='F').copy()
            elif '_y' in field:
                data = data[1::3].reshape(grid_dims,order='F').copy()
            elif '_z' in field:
                data = data[2::3].reshape(grid_dims,order='F').copy()
        f.close()
        if grid.pf.field_ordering == 1:
            return data.T.astype("float64")
        else:
            return data.astype("float64")

    def _read_data_slice(self, grid, field, axis, coord):
        sl = [slice(None), slice(None), slice(None)]
        sl[axis] = slice(coord, coord + 1)
        if grid.pf.field_ordering == 1:
            sl.reverse()
        return self._read_data_set(grid, field)[sl]


def get_read_table_offset(f):
    line = f.readline()
    while True:
        splitup = line.strip().split()
        if 'CELL_DATA' in splitup:
            f.readline()
            read_table_offset = f.tell()
            break
        line = f.readline()
    return read_table_offset


