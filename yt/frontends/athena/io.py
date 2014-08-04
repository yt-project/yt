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
from yt.funcs import mylog, defaultdict

class IOHandlerAthena(BaseIOHandler):
    _dataset_type = "athena"
    _offset_string = 'data:offsets=0'
    _data_string = 'data:datatype=0'
    _read_table_offset = None

    def _field_dict(self,fhandle):
        keys = fhandle['field_types'].keys()
        val = fhandle['field_types'].keys()
        return dict(zip(keys,val))

    def _read_field_names(self,grid):
        pass

    def _read_chunk_data(self,chunk,fields):
        data = {}
        grids_by_file = defaultdict(list)
        if len(chunk.objs) == 0: return data
        field_list = set(f[1] for f in fields)
        for grid in chunk.objs:
            if grid.filename is None:
                continue
            f = open(grid.filename, "rb")
            data[grid.id] = {}
            grid_ncells = np.prod(grid.ActiveDimensions)
            grid_dims = grid.ActiveDimensions
            grid0_ncells = np.prod(grid.index.grid_dimensions[0,:])
            read_table_offset = get_read_table_offset(f)
            for field in self.ds.field_list:
                dtype, offsetr = grid.index._field_map[field]
                if grid_ncells != grid0_ncells:
                    offset = offsetr + ((grid_ncells-grid0_ncells) * (offsetr//grid0_ncells))
                if grid_ncells == grid0_ncells:
                    offset = offsetr
                f.seek(read_table_offset+offset)
                if dtype == 'scalar':
                    v = np.fromfile(f, dtype='>f4',
                                    count=grid_ncells).reshape(grid_dims,order='F')
                if dtype == 'vector':
                    v = np.fromfile(f, dtype='>f4', count=3*grid_ncells)
                if '_x' in field[-1]:
                    v = v[0::3].reshape(grid_dims,order='F')
                elif '_y' in field[-1]:
                    v = v[1::3].reshape(grid_dims,order='F')
                elif '_z' in field[-1]:
                    v = v[2::3].reshape(grid_dims,order='F')
                if grid.ds.field_ordering == 1:
                    data[grid.id][field] = v.T.astype("float64")
                else:
                    data[grid.id][field] = v.astype("float64")
            f.close()
        return data
    
    def _read_data_slice(self, grid, field, axis, coord):
        sl = [slice(None), slice(None), slice(None)]
        sl[axis] = slice(coord, coord + 1)
        if grid.ds.field_ordering == 1:
            sl.reverse()
        return self._read_data_set(grid, field)[sl]

    def _read_fluid_selection(self, chunks, selector, fields, size):
        chunks = list(chunks)
        if any((ftype != "athena" for ftype, fname in fields)):
            raise NotImplementedError
        rv = {}
        for field in fields:
            rv[field] = np.empty(size, dtype="float64")
        ng = sum(len(c.objs) for c in chunks)
        mylog.debug("Reading %s cells of %s fields in %s grids",
                    size, [f2 for f1, f2 in fields], ng)
        ind = 0
        for chunk in chunks:
            data = self._read_chunk_data(chunk, fields)
            for g in chunk.objs:
                for field in fields:
                    ftype, fname = field
                    ds = data[g.id].pop(field)
                    nd = g.select(selector, ds, rv[field], ind) # caches
                ind += nd
                data.pop(g.id)
        return rv

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


