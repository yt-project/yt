"""
GAMER-specific IO functions



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2016, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
from yt.funcs import mylog
from yt.utilities.io_handler import \
    BaseIOHandler

class IOHandlerGAMER(BaseIOHandler):
    _particle_reader = False
    _dataset_type    = "gamer"

    def __init__(self, ds):
        super(IOHandlerGAMER, self).__init__(ds)
        self._handle      = ds._handle
        self._field_dtype = "float64" # fixed even when FLOAT8 is off

    def _read_particle_coords(self, chunks, ptf):
        # This needs to *yield* a series of tuples of (ptype, (x, y, z)).
        # chunks is a list of chunks, and ptf is a dict where the keys are
        # ptypes and the values are lists of fields.
        pass

    def _read_particle_fields(self, chunks, ptf, selector):
        # This gets called after the arrays have been allocated.  It needs to
        # yield ((ptype, field), data) where data is the masked results of
        # reading ptype, field and applying the selector to the data read in.
        # Selector objects have a .select_points(x,y,z) that returns a mask, so
        # you need to do your masking here.
        pass


    def _read_fluid_selection(self, chunks, selector, fields, size):
        chunks = list(chunks) # generator --> list

        if any( (ftype != "gamer" for ftype, fname in fields) ):
            raise NotImplementedError

        rv = {}
        for field in fields: rv[field] = np.empty( size, dtype=self._field_dtype )

        ng = sum( len(c.objs) for c in chunks ) # c.objs is a list of grids
        mylog.debug( "Reading %s cells of %s fields in %s grids",
                     size, [f2 for f1, f2 in fields], ng )
        offset = 0
        for chunk in chunks:
            data = self._read_chunk_data( chunk, fields )
            for g in chunk.objs:
                for field in fields:
                    ds    = data[g.id].pop(field)
                    # array return from g.select (i.e., rv[field]) is flat
                    ncell = g.select( selector, ds, rv[field], offset )
                offset += ncell
                data.pop(g.id)
        return rv

    def _read_chunk_data(self, chunk, fields):
        data = {}
        if len(chunk.objs) == 0: return data

        for g in chunk.objs:
            data[g.id] = {}
            LvName = 'Level_%02i'   % g.Level
            CName  = 'Cluster_%09i' % g.CID
            PName  = 'Patch_%09i'   % g.PID

            for field in fields:
                # transpose x-z since YT assumes that consecutive cells along z
                # are contiguous in memory
                data[g.id][field] \
                = self._handle[LvName][CName][PName][ field[1] ].swapaxes(0,2)
        return data

