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
from yt.geometry.selection_routines import AlwaysSelector


#-----------------------------------------------------------------------------
# GAMER shares a similar HDF5 format, and thus io.py as well, with FLASH
#-----------------------------------------------------------------------------


# group grids with consecutive indices together to improve the I/O performance
# --> grids are assumed to be sorted into ascending numerical order already
def grid_sequences(grids):
    for k, g in groupby( enumerate(grids), lambda i_x:i_x[0]-i_x[1].id ):
        seq = list(v[1] for v in g)
        yield seq

def particle_sequences(grids):
    for k, g in groupby( enumerate(grids), lambda i_x:i_x[0]-i_x[1].id ):
        seq = list(v[1] for v in g)
        yield seq[0], seq[-1]

class IOHandlerGAMER(BaseIOHandler):
    _particle_reader = False
    _dataset_type    = "gamer"

    def __init__(self, ds):
        super(IOHandlerGAMER, self).__init__(ds)
        self._handle          = ds._handle
        self._field_dtype     = "float64" # fixed even when FLOAT8 is off
        self._particle_handle = ds._particle_handle

    def _read_particle_coords(self, chunks, ptf):
        chunks = list(chunks)   # generator --> list
        p_idx  = self.ds.index._particle_indices

        # shortcuts
        par_posx = self._handle[ "/Particle/ParPosX" ]
        par_posy = self._handle[ "/Particle/ParPosY" ]
        par_posz = self._handle[ "/Particle/ParPosZ" ]

        # currently GAMER does not support multiple particle types
        assert( len(ptf) == 1 )
        ptype = ptf.keys()[0]

        for chunk in chunks:
            for g1, g2 in particle_sequences(chunk.objs):
                start = p_idx[g1.id    ]
                end   = p_idx[g2.id + 1]
                x     = np.asarray( par_posx[start:end], dtype=self._field_dtype )
                y     = np.asarray( par_posy[start:end], dtype=self._field_dtype )
                z     = np.asarray( par_posz[start:end], dtype=self._field_dtype )
                yield ptype, (x, y, z)

    def _read_particle_fields(self, chunks, ptf, selector):
        chunks = list(chunks)   # generator --> list
        p_idx  = self.ds.index._particle_indices

        # shortcuts
        par_posx = self._handle[ "/Particle/ParPosX" ]
        par_posy = self._handle[ "/Particle/ParPosY" ]
        par_posz = self._handle[ "/Particle/ParPosZ" ]

        # currently GAMER does not support multiple particle types
        assert( len(ptf) == 1 )
        ptype   = ptf.keys()[0]
        pfields = ptf[ptype]

        for chunk in chunks:
            for g1, g2 in particle_sequences(chunk.objs):
                start = p_idx[g1.id    ]
                end   = p_idx[g2.id + 1]
                x     = np.asarray( par_posx[start:end], dtype=self._field_dtype )
                y     = np.asarray( par_posy[start:end], dtype=self._field_dtype )
                z     = np.asarray( par_posz[start:end], dtype=self._field_dtype )

                mask = selector.select_points(x, y, z, 0.0)
                if mask is None: continue

                for field in pfields:
                    data = self._handle["/Particle/%s" % field][start:end]
                    yield (ptype, field), data[mask]

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
            try:
                ds = self._handle[ "/GridData/%s" % field[1] ]
            except KeyError:
                ds = self._handle[ "/Data/%s" % field[1] ]

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

        # Split into particles and non-particles
        fluid_fields, particle_fields = [], []
        for ftype, fname in fields:
            if ftype in self.ds.particle_types:
                particle_fields.append( (ftype, fname) )
            else:
                fluid_fields.append( (ftype, fname) )

        # particles
        if len(particle_fields) > 0:
            selector = AlwaysSelector(self.ds)
            rv.update( self._read_particle_selection(
                [chunk], selector, particle_fields) )

        # fluid
        if len(fluid_fields) == 0: return rv

        for field in fluid_fields:
            try:
                ds = self._handle[ "/GridData/%s" % field[1] ]
            except KeyError:
                ds = self._handle[ "/Data/%s" % field[1] ]

            for gs in grid_sequences(chunk.objs):
                start = gs[ 0].id
                end   = gs[-1].id + 1
                data  = ds[start:end,:,:,:].transpose()
                for i, g in enumerate(gs):
                    rv[g.id][field] = np.asarray( data[...,i], dtype=self._field_dtype )
        return rv
