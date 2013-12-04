"""
FLASH-specific IO functions



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
import h5py
from itertools import groupby

from yt.utilities.io_handler import \
    BaseIOHandler
from yt.utilities.logger import ytLogger as mylog

# http://stackoverflow.com/questions/2361945/detecting-consecutive-integers-in-a-list
def particle_sequences(grids):
    g_iter = sorted(grids, key = lambda g: g.id)
    for k, g in groupby(enumerate(g_iter), lambda (i,x):i-x.id):
        seq = list(v[1] for v in g)
        yield seq[0], seq[-1]

def grid_sequences(grids):
    g_iter = sorted(grids, key = lambda g: g.id)
    for k, g in groupby(enumerate(g_iter), lambda (i,x):i-x.id):
        seq = list(v[1] for v in g)
        yield seq

class IOHandlerFLASH(BaseIOHandler):
    _particle_reader = False
    _data_style = "flash_hdf5"

    def __init__(self, pf):
        super(IOHandlerFLASH, self).__init__(pf)
        # Now we cache the particle fields
        self._handle = pf._handle
        self._particle_handle = pf._particle_handle
        
        try :
            particle_fields = [s[0].strip() for s in
                               self._particle_handle["/particle names"][:]]
            self._particle_fields = dict([("particle_" + s, i) for i, s in
                                          enumerate(particle_fields)])
        except KeyError:
            self._particle_fields = {}

    def _read_particles(self, fields_to_read, type, args, grid_list,
            count_list, conv_factors):
        pass

    def _read_particle_coords(self, chunks, ptf):
        chunks = list(chunks)
        f_part = self._particle_handle
        p_ind = self.pf.h._particle_indices
        px, py, pz = (self._particle_fields["particle_pos%s" % ax]
                      for ax in 'xyz')
        p_fields = f_part["/tracer particles"]
        assert(len(ptf) == 1)
        ptype = ptf.keys()[0]
        for chunk in chunks:
            start = end = None
            for g1, g2 in particle_sequences(chunk.objs):
                start = p_ind[g1.id - g1._id_offset]
                end = p_ind[g2.id - g2._id_offset + 1]
                x = p_fields[start:end, px]
                y = p_fields[start:end, py]
                z = p_fields[start:end, pz]
                yield ptype, (x, y, z)

    def _read_particle_fields(self, chunks, ptf, selector):
        chunks = list(chunks)
        f_part = self._particle_handle
        p_ind = self.pf.h._particle_indices
        px, py, pz = (self._particle_fields["particle_pos%s" % ax]
                      for ax in 'xyz')
        p_fields = f_part["/tracer particles"]
        assert(len(ptf) == 1)
        ptype = ptf.keys()[0]
        field_list = ptf[ptype]
        for chunk in chunks:
            for g1, g2 in particle_sequences(chunk.objs):
                start = p_ind[g1.id - g1._id_offset]
                end = p_ind[g2.id - g2._id_offset + 1]
                x = p_fields[start:end, px]
                y = p_fields[start:end, py]
                z = p_fields[start:end, pz]
                mask = selector.select_points(x, y, z)
                if mask is None: continue
                for field in field_list:
                    fi = self._particle_fields[field]
                    data = p_fields[start:end, fi]
                    yield (ptype, field), data[mask]

    def _read_data_set(self, grid, field):
        f = self._handle
        f_part = self._particle_handle
        if field in self._particle_fields:
            if grid.NumberOfParticles == 0: return np.array([], dtype='float64')
            start = self.pf.h._particle_indices[grid.id - grid._id_offset]
            end = self.pf.h._particle_indices[grid.id - grid._id_offset + 1]
            fi = self._particle_fields[field]
            tr = f_part["/tracer particles"][start:end, fi]
        else:
            tr = f["/%s" % field][grid.id - grid._id_offset,:,:,:].transpose()
        return tr.astype("float64")

    def _read_data_slice(self, grid, field, axis, coord):
        sl = [slice(None), slice(None), slice(None)]
        sl[axis] = slice(coord, coord + 1)
        f = self._handle
        tr = f["/%s" % field][grid.id - grid._id_offset].transpose()[sl]
        return tr.astype("float64")

    def _read_fluid_selection(self, chunks, selector, fields, size):
        chunks = list(chunks)
        if any((ftype != "gas" for ftype, fname in fields)):
            raise NotImplementedError
        f = self._handle
        rv = {}
        for field in fields:
            ftype, fname = field
            dt = f["/%s" % fname].dtype
            if dt == "float32": dt = "float64"
            rv[field] = np.empty(size, dtype=dt)
        ng = sum(len(c.objs) for c in chunks)
        mylog.debug("Reading %s cells of %s fields in %s blocks",
                    size, [f2 for f1, f2 in fields], ng)
        for field in fields:
            ftype, fname = field
            ds = f["/%s" % fname]
            ind = 0
            for chunk in chunks:
                for gs in grid_sequences(chunk.objs):
                    start = gs[0].id - gs[0]._id_offset
                    end = gs[-1].id - gs[-1]._id_offset + 1
                    data = ds[start:end,:,:,:].transpose()
                    for i, g in enumerate(gs):
                        ind += g.select(selector, data[...,i], rv[field], ind)
        return rv

    def _read_chunk_data(self, chunk, fields):
        f = self._handle
        rv = {}
        for g in chunk.objs:
            rv[g.id] = {}
        for field in fields:
            ftype, fname = field
            ds = f["/%s" % fname]
            ind = 0
            for gs in grid_sequences(chunk.objs):
                start = gs[0].id - gs[0]._id_offset
                end = gs[-1].id - gs[-1]._id_offset + 1
                data = ds[start:end,:,:,:].transpose()
                for i, g in enumerate(gs):
                    rv[g.id][field] = data[...,i]
        return rv

