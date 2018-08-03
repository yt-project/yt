"""
The data-file handling functions



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2018, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.utilities.io_handler import \
    BaseIOHandler
import numpy as np
from yt.funcs import mylog

class IOHandlerCholla(BaseIOHandler):
    _particle_reader = False
    _dataset_type = "cholla"

    def __init__(self,ds):
        super(IOHandlerCholla, self).__init__(ds)
        self._handle = ds._handle

    def _read_particles(self, fields_to_read, type, args, grid_list, 
                        count_list, conv_factors):
        pass

    def _read_fluid_selection(self, chunks, selector, fields, size):
        chunks = list(chunks)
        if any((ftype != "cholla" for ftype, fname in fields)):
            raise NotImplementedError
        f = self._handle
        rv = {}
        for field in fields:
            # Always use *native* 64-bit float.
            rv[field] = np.empty(size, dtype="=f8")
        ng = sum(len(c.objs) for c in chunks)
        mylog.debug("Reading %s cells of %s fields in %s blocks",
                    size, [f2 for f1, f2 in fields], ng)
        for field in fields:
            ftype, fname = field
            dname, fdi = self.ds._field_map[fname]
            ds = f["/%s" % dname]
            ind = 0
            for chunk in chunks:
                if self.ds.logarithmic:
                    for mesh in chunk.objs:
                        nx, ny, nz = mesh.mesh_dims // self.ds.index.mesh_factors
                        data = np.empty(mesh.mesh_dims, dtype="=f8")
                        for n, id in enumerate(mesh.mesh_blocks):
                            data[ii[n]*nx:(ii[n]+1)*nx,jj[n]*ny:(jj[n]+1)*ny,kk[n]*nz:(kk[n]+1)*nz] = \
                                 ds[fdi,id,:,:,:].transpose()
                        ind += mesh.select(selector, data, rv[field], ind)  # caches
                else:
                    for gs in grid_sequences(chunk.objs):
                        start = gs[0].id - gs[0]._id_offset
                        end = gs[-1].id - gs[-1]._id_offset + 1
                        data = ds[fdi,start:end,:,:,:].transpose()
                        for i, g in enumerate(gs):
                            ind += g.select(selector, data[...,i], rv[field], ind)
        return rv
