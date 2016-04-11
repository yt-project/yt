"""
Athena++-specific IO functions



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
from itertools import groupby

from yt.utilities.io_handler import \
    BaseIOHandler
from yt.utilities.logger import ytLogger as mylog

# http://stackoverflow.com/questions/2361945/detecting-consecutive-integers-in-a-list
def grid_sequences(grids):
    g_iter = sorted(grids, key = lambda g: g.id)
    for k, g in groupby(enumerate(g_iter), lambda i_x1:i_x1[0]-i_x1[1].id):
        seq = list(v[1] for v in g)
        yield seq

class IOHandlerAthenaPP(BaseIOHandler):
    _particle_reader = False
    _dataset_type = "athena++"

    def __init__(self, ds):
        super(IOHandlerAthenaPP, self).__init__(ds)
        self._handle = ds._handle

    def _read_particles(self, fields_to_read, type, args, grid_list,
            count_list, conv_factors):
        pass

    def _read_fluid_selection(self, chunks, selector, fields, size):
        chunks = list(chunks)
        if any((ftype != "athena++" for ftype, fname in fields)):
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
            f_fname = self.ds._field_map[field[1]]
            ind = 0
            for chunk in chunks:
                for gs in grid_sequences(chunk.objs):
                    for i, g in enumerate(gs):
                        data = f["MeshBlock%d" % g.id][f_fname][:,:,:].transpose()
                        ind += g.select(selector, data, rv[field], ind)
        return rv

    def _read_chunk_data(self, chunk, fields):
        f = self._handle
        rv = {}
        for g in chunk.objs:
            rv[g.id] = {}
        if len(fields) == 0:
            return rv
        for field in fields:
            f_fname = self.ds._field_map[field[1]]
            for gs in grid_sequences(chunk.objs):
                for i, g in enumerate(gs):
                    data = f["MeshBlock%d" % g.id][f_fname][:,:,:].transpose()
                    rv[g.id][field] = np.asarray(data, "=f8")
        return rv

