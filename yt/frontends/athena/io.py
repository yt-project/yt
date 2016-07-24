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
from yt.funcs import mylog
from .data_structures import chk23

float_size = {"float":np.dtype(">f4").itemsize,
              "double":np.dtype(">f8").itemsize}

axis_list = ["_x","_y","_z"]

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
        if len(chunk.objs) == 0: return data
        for grid in chunk.objs:
            if grid.filename is None:
                continue
            f = open(grid.filename, "rb")
            data[grid.id] = {}
            grid_dims = grid.ActiveDimensions
            read_dims = grid.read_dims.astype("int64")
            grid_ncells = np.prod(read_dims)
            grid0_ncells = np.prod(grid.index.grids[0].read_dims)
            read_table_offset = get_read_table_offset(f)
            for field in fields:
                ftype, offsetr, dtype = grid.index._field_map[field]
                if grid_ncells != grid0_ncells:
                    offset = offsetr + ((grid_ncells-grid0_ncells) * (offsetr//grid0_ncells))
                if grid_ncells == grid0_ncells:
                    offset = offsetr
                offset = int(offset) # Casting to be certain.
                file_offset = grid.file_offset[2]*read_dims[0]*read_dims[1]*float_size[dtype]
                xread = slice(grid.file_offset[0],grid.file_offset[0]+grid_dims[0])
                yread = slice(grid.file_offset[1],grid.file_offset[1]+grid_dims[1])
                f.seek(read_table_offset+offset+file_offset)
                if dtype == 'float':
                    dt = '>f4'
                elif dtype == 'double':
                    dt = '>f8'
                if ftype == 'scalar':
                    f.seek(read_table_offset+offset+file_offset)
                    v = np.fromfile(f, dtype=dt,
                                    count=grid_ncells).reshape(read_dims,order='F')
                if ftype == 'vector':
                    vec_offset = axis_list.index(field[-1][-2:])
                    f.seek(read_table_offset+offset+3*file_offset)
                    v = np.fromfile(f, dtype=dt, count=3*grid_ncells)
                    v = v[vec_offset::3].reshape(read_dims,order='F')
                if grid.ds.field_ordering == 1:
                    data[grid.id][field] = v[xread,yread,:].T.astype("float64")
                else:
                    data[grid.id][field] = v[xread,yread,:].astype("float64")
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
        chkc = chk23('CELL_DATA')
        chkp = chk23('POINT_DATA')
        if chkc in splitup or chkp in splitup:
            f.readline()
            read_table_offset = f.tell()
            break
        line = f.readline()
    return read_table_offset


