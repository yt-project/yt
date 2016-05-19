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
from itertools import groupby

from yt.utilities.io_handler import \
    BaseIOHandler
from yt.utilities.logger import ytLogger as mylog
from yt.geometry.selection_routines import AlwaysSelector
from yt.utilities.lib.geometry_utils import \
    compute_morton

# http://stackoverflow.com/questions/2361945/detecting-consecutive-integers-in-a-list
def particle_sequences(grids):
    g_iter = sorted(grids, key = lambda g: g.id)
    for k, g in groupby(enumerate(g_iter), lambda i_x:i_x[0]-i_x[1].id):
        seq = list(v[1] for v in g)
        yield seq[0], seq[-1]

def grid_sequences(grids):
    g_iter = sorted(grids, key = lambda g: g.id)
    for k, g in groupby(enumerate(g_iter), lambda i_x1:i_x1[0]-i_x1[1].id):
        seq = list(v[1] for v in g)
        yield seq

def determine_particle_fields(handle):
    try:
        particle_fields = [s[0].decode("ascii","ignore").strip()
                           for s in handle["/particle names"][:]]
        _particle_fields = dict([("particle_" + s, i) for i, s in
                                 enumerate(particle_fields)])
    except KeyError:
        _particle_fields = {}
    return _particle_fields

class IOHandlerFLASH(BaseIOHandler):
    _particle_reader = False
    _dataset_type = "flash_hdf5"

    def __init__(self, ds):
        super(IOHandlerFLASH, self).__init__(ds)
        # Now we cache the particle fields
        self._handle = ds._handle
        self._particle_handle = ds._particle_handle
        self._particle_fields = determine_particle_fields(self._particle_handle)

    def _read_particles(self, fields_to_read, type, args, grid_list,
            count_list, conv_factors):
        pass

    def _read_particle_coords(self, chunks, ptf):
        chunks = list(chunks)
        f_part = self._particle_handle
        p_ind = self.ds.index._particle_indices
        px, py, pz = (self._particle_fields["particle_pos%s" % ax]
                      for ax in 'xyz')
        p_fields = f_part["/tracer particles"]
        assert(len(ptf) == 1)
        ptype = list(ptf.keys())[0]
        for chunk in chunks:
            start = end = None
            for g1, g2 in particle_sequences(chunk.objs):
                start = p_ind[g1.id - g1._id_offset]
                end = p_ind[g2.id - g2._id_offset + 1]
                x = np.asarray(p_fields[start:end, px], dtype="=f8")
                y = np.asarray(p_fields[start:end, py], dtype="=f8")
                z = np.asarray(p_fields[start:end, pz], dtype="=f8")
                yield ptype, (x, y, z)

    def _read_particle_fields(self, chunks, ptf, selector):
        chunks = list(chunks)
        f_part = self._particle_handle
        p_ind = self.ds.index._particle_indices
        px, py, pz = (self._particle_fields["particle_pos%s" % ax]
                      for ax in 'xyz')
        p_fields = f_part["/tracer particles"]
        assert(len(ptf) == 1)
        ptype = list(ptf.keys())[0]
        field_list = ptf[ptype]
        for chunk in chunks:
            for g1, g2 in particle_sequences(chunk.objs):
                start = p_ind[g1.id - g1._id_offset]
                end = p_ind[g2.id - g2._id_offset + 1]
                x = np.asarray(p_fields[start:end, px], dtype="=f8")
                y = np.asarray(p_fields[start:end, py], dtype="=f8")
                z = np.asarray(p_fields[start:end, pz], dtype="=f8")
                mask = selector.select_points(x, y, z, 0.0)
                if mask is None: continue
                for field in field_list:
                    fi = self._particle_fields[field]
                    data = p_fields[start:end, fi]
                    yield (ptype, field), data[mask]

    def _read_fluid_selection(self, chunks, selector, fields, size):
        chunks = list(chunks)
        if any((ftype != "flash" for ftype, fname in fields)):
            raise NotImplementedError
        f = self._handle
        rv = {}
        for field in fields:
            ftype, fname = field
            # Always use *native* 64-bit float.
            rv[field] = np.empty(size, dtype="=f8")
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
        # Split into particles and non-particles
        fluid_fields, particle_fields = [], []
        for ftype, fname in fields:
            if ftype in self.ds.particle_types:
                particle_fields.append((ftype, fname))
            else:
                fluid_fields.append((ftype, fname))
        if len(particle_fields) > 0:
            selector = AlwaysSelector(self.ds)
            rv.update(self._read_particle_selection(
                [chunk], selector, particle_fields))
        if len(fluid_fields) == 0: return rv
        for field in fluid_fields:
            ftype, fname = field
            ds = f["/%s" % fname]
            for gs in grid_sequences(chunk.objs):
                start = gs[0].id - gs[0]._id_offset
                end = gs[-1].id - gs[-1]._id_offset + 1
                data = ds[start:end,:,:,:].transpose()
                for i, g in enumerate(gs):
                    rv[g.id][field] = np.asarray(data[...,i], "=f8")
        return rv

class IOHandlerFLASHParticle(BaseIOHandler):
    _particle_reader = True
    _dataset_type = "flash_particle_hdf5"

    def __init__(self, ds):
        super(IOHandlerFLASHParticle, self).__init__(ds)
        # Now we cache the particle fields
        self._handle = ds._handle
        self._particle_fields = determine_particle_fields(self._handle)
        self._position_fields = [self._particle_fields["particle_pos%s" % ax]
                                 for ax in 'xyz']
        self._chunksize = 32**3

    def _read_fluid_selection(self, chunks, selector, fields, size):
        raise NotImplementedError

    def _read_particle_coords(self, chunks, ptf):
        chunks = list(chunks)
        data_files = set([])
        assert(len(ptf) == 1)
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        px, py, pz = self._position_fields
        p_fields = self._handle["/tracer particles"]
        assert(len(data_files) == 1)
        for data_file in sorted(data_files):
            pcount = self._count_particles(data_file)["io"]
            for ptype, field_list in sorted(ptf.items()):
                total = 0
                while total < pcount:
                    count = min(self._chunksize, pcount - total)
                    x = np.asarray(p_fields[total:total+count, px], dtype="=f8")
                    y = np.asarray(p_fields[total:total+count, py], dtype="=f8")
                    z = np.asarray(p_fields[total:total+count, pz], dtype="=f8")
                    total += count
                    yield ptype, (x, y, z)

    def _read_particle_fields(self, chunks, ptf, selector):
        chunks = list(chunks)
        data_files = set([])
        assert(len(ptf) == 1)
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        px, py, pz = self._position_fields
        p_fields = self._handle["/tracer particles"]
        assert(len(data_files) == 1)
        for data_file in sorted(data_files):
            pcount = self._count_particles(data_file)["io"]
            for ptype, field_list in sorted(ptf.items()):
                total = 0
                while total < pcount:
                    count = min(self._chunksize, pcount - total)
                    x = np.asarray(p_fields[total:total+count, px], dtype="=f8")
                    y = np.asarray(p_fields[total:total+count, py], dtype="=f8")
                    z = np.asarray(p_fields[total:total+count, pz], dtype="=f8")
                    total += count
                    mask = selector.select_points(x, y, z, 0.0)
                    del x, y, z
                    if mask is None: continue
                    for field in field_list:
                        fi = self._particle_fields[field]
                        data = p_fields[total-count:total, fi]
                        yield (ptype, field), data[mask]

    def _initialize_index(self, data_file, regions):
        p_fields = self._handle["/tracer particles"]
        px, py, pz = self._position_fields
        pcount = self._count_particles(data_file)["io"]
        morton = np.empty(pcount, dtype='uint64')
        ind = 0
        while ind < pcount:
            npart = min(self._chunksize, pcount - ind)
            pos = np.empty((npart, 3), dtype="=f8")
            pos[:,0] = p_fields[ind:ind+npart, px]
            pos[:,1] = p_fields[ind:ind+npart, py]
            pos[:,2] = p_fields[ind:ind+npart, pz]
            regions.add_data_file(pos, data_file.file_id)
            morton[ind:ind+npart] = \
                compute_morton(pos[:,0], pos[:,1], pos[:,2],
                               data_file.ds.domain_left_edge,
                               data_file.ds.domain_right_edge)
            ind += self._chunksize
        return morton

    def _count_particles(self, data_file):
        pcount = {"io": self._handle["/localnp"][:].sum()}
        return pcount

    def _identify_fields(self, data_file):
        fields = [("io", field) for field in self._particle_fields]
        return fields, {}
