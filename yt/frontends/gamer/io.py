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
from itertools import groupby

from yt.utilities.io_handler import \
    BaseIOHandler
from yt.utilities.logger import ytLogger as mylog


#-----------------------------------------------------------------------------
# GAMER shares a similar HDF5 format, and thus io.py as well, with FLASH
#-----------------------------------------------------------------------------


# group grids with consecutive indices together to improve the I/O performance
def grid_sequences(grids):
    for k, g in groupby( enumerate(grids), lambda i_x1:i_x1[0]-i_x1[1].id ):
        seq = list(v[1] for v in g)
        yield seq

class IOHandlerGAMER(BaseIOHandler):
    _particle_reader = False
    _dataset_type    = "gamer"

    def __init__(self, ds):
        super(IOHandlerGAMER, self).__init__(ds)
        self._handle      = ds._handle
        self._field_dtype = "float64" # fixed even when FLOAT8 is off

    def _read_particle_coords(self, chunks, ptf):
        pass

    def _read_particle_fields(self, chunks, ptf, selector):
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

        for field in fields:
            ds     = self._handle[ "/Data/%s" % field[1] ]
            offset = 0
            for chunk in chunks:
                for gs in grid_sequences(chunk.objs):
                    start = gs[ 0].id
                    end   = gs[-1].id + 1
                    data  = ds[start:end,:,:,:].transpose()
                    for i, g in enumerate(gs):
                        offset += g.select( selector, data[...,i], rv[field], offset )
        return rv

    def _read_chunk_data(self, chunk, fields):
        rv = {}
        if len(chunk.objs) == 0: return rv 

        for g in chunk.objs: rv[g.id] = {}

        for field in fields:
            ds = self._handle[ "/Data/%s" % field[1] ]

            for gs in grid_sequences(chunk.objs):
                start = gs[ 0].id
                end   = gs[-1].id + 1
                data  = ds[start:end,:,:,:].transpose()
                for i, g in enumerate(gs):
                    rv[g.id][field] = np.asarray( data[...,i], dtype=self._field_dtype )
        return rv
