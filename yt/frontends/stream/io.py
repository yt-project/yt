"""
Enzo-specific IO functions



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from yt.utilities.io_handler import \
    BaseIOHandler
from yt.utilities.logger import ytLogger as mylog
from yt.utilities.lib.geometry_utils import compute_morton
from yt.utilities.exceptions import YTDomainOverflow

class IOHandlerStream(BaseIOHandler):

    _dataset_type = "stream"
    _vector_fields = ("particle_velocity", "particle_position")

    def __init__(self, ds):
        self.fields = ds.stream_handler.fields
        self.field_units = ds.stream_handler.field_units
        super(IOHandlerStream, self).__init__(ds)

    def _read_data_set(self, grid, field):
        # This is where we implement processor-locking
        #if grid.id not in self.grids_in_memory:
        #    mylog.error("Was asked for %s but I have %s", grid.id, self.grids_in_memory.keys())
        #    raise KeyError
        tr = self.fields[grid.id][field]
        # If it's particles, we copy.
        if len(tr.shape) == 1: return tr.copy()
        # New in-place unit conversion breaks if we don't copy first
        return tr

    def _read_fluid_selection(self, chunks, selector, fields, size):
        chunks = list(chunks)
        if any((ftype not in self.ds.fluid_types for ftype, fname in fields)):
            raise NotImplementedError
        rv = {}
        for field in fields:
            rv[field] = self.ds.arr(np.empty(size, dtype="float64"))
        ng = sum(len(c.objs) for c in chunks)
        mylog.debug("Reading %s cells of %s fields in %s blocks",
                    size, [f2 for f1, f2 in fields], ng)
        for field in fields:
            ftype, fname = field
            ind = 0
            for chunk in chunks:
                for g in chunk.objs:
                    ds = self.fields[g.id][ftype, fname]
                    ind += g.select(selector, ds, rv[field], ind) # caches
        return rv

    def _read_particle_coords(self, chunks, ptf):
        chunks = list(chunks)
        for chunk in chunks:
            for g in chunk.objs:
                if g.NumberOfParticles == 0: continue
                gf = self.fields[g.id]
                for ptype, field_list in sorted(ptf.items()):
                    if (ptype, "particle_position") in gf:
                        x, y, z = gf[ptype, "particle_position"].T
                    else:
                        x, y, z = (gf[ptype, "particle_position_%s" % ax] for
                                   ax in 'xyz')
                    yield ptype, (x, y, z)

    def _read_particle_fields(self, chunks, ptf, selector):
        chunks = list(chunks)
        for chunk in chunks:
            for g in chunk.objs:
                if g.NumberOfParticles == 0: continue
                gf = self.fields[g.id]
                for ptype, field_list in sorted(ptf.items()):
                    if (ptype, "particle_position") in gf:
                        x, y, z = gf[ptype, "particle_position"].T
                    else:
                        x, y, z = (gf[ptype, "particle_position_%s" % ax] for
                                   ax in 'xyz')
                    mask = selector.select_points(x, y, z, 0.0)
                    if mask is None: continue
                    for field in field_list:
                        data = np.asarray(gf[ptype, field])
                        yield (ptype, field), data[mask]

    @property
    def _read_exception(self):
        return KeyError

class StreamParticleIOHandler(BaseIOHandler):

    _vector_fields = ("particle_position", "particle_velocity")
    _dataset_type = "stream_particles"
    _vector_fields = ("particle_velocity", "particle_position")

    def __init__(self, ds):
        self.fields = ds.stream_handler.fields
        super(StreamParticleIOHandler, self).__init__(ds)

    def _read_particle_coords(self, chunks, ptf):
        chunks = list(chunks)
        data_files = set([])
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        for data_file in sorted(data_files):
            f = self.fields[data_file.filename]
            # This double-reads
            for ptype, field_list in sorted(ptf.items()):
                yield ptype, (f[ptype, "particle_position_x"],
                              f[ptype, "particle_position_y"],
                              f[ptype, "particle_position_z"])
            
    def __count_particles_chunks(self, chunks, ptf, selector):
        # DISABLED
        # I have left this in here, but disabled, because of two competing
        # problems:
        #   * The IndexedOctreeSubsetSelector currently selects *all* particles
        #   * Slicing a deposited field thus throws an error, since the octree
        #     estimate fails.
        #   * BUT, it provides considerable speedup in some situations for
        #     stream datasets.
        # So, pending its re-enabling, we'll leave it here.
        # 
        # This is allowed to over-estimate.  We probably *will*, too, because
        # we're going to count *all* of the particles, not just individual
        # types.
        count = 0
        psize = {}
        for chunk in chunks:
            for obj in chunk.objs:
                count += selector.count_octs(obj.oct_handler, obj.domain_id)
        for ptype in ptf:
            psize[ptype] = self.ds.n_ref * count
        return psize

    def _read_particle_fields(self, chunks, ptf, selector):
        data_files = set([])
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        for data_file in sorted(data_files):
            f = self.fields[data_file.filename]
            for ptype, field_list in sorted(ptf.items()):
                if (ptype, "particle_position") in f:
                    x = f[ptype, "particle_position"][:,0]
                    y = f[ptype, "particle_position"][:,1]
                    z = f[ptype, "particle_position"][:,2]
                else:
                    x, y, z = (f[ptype, "particle_position_%s" % ax]
                               for ax in 'xyz')
                mask = selector.select_points(x, y, z, 0.0)
                if mask is None: continue
                for field in field_list:
                    data = f[ptype, field][mask]
                    yield (ptype, field), data

    def _initialize_index(self, data_file, regions):
        # self.fields[g.id][fname] is the pattern here
        morton = []
        for ptype in self.ds.particle_types_raw:
            try:
                pos = np.column_stack(self.fields[data_file.filename][
                    (ptype, "particle_position_%s" % ax)] for ax in 'xyz')
            except KeyError:
                pos = self.fields[data_file.filename][ptype, "particle_position"]
            if np.any(pos.min(axis=0) < data_file.ds.domain_left_edge) or \
               np.any(pos.max(axis=0) > data_file.ds.domain_right_edge):
                raise YTDomainOverflow(pos.min(axis=0), pos.max(axis=0),
                                       data_file.ds.domain_left_edge,
                                       data_file.ds.domain_right_edge)
            regions.add_data_file(pos, data_file.file_id)
            morton.append(compute_morton(
                    pos[:,0], pos[:,1], pos[:,2],
                    data_file.ds.domain_left_edge,
                    data_file.ds.domain_right_edge))
        return np.concatenate(morton)

    def _count_particles(self, data_file):
        pcount = {}
        for ptype in self.ds.particle_types_raw:
            d = self.fields[data_file.filename]
            try:
                pcount[ptype] = d[ptype, "particle_position_x"].size
            except KeyError:
                pcount[ptype] = d[ptype, "particle_position"].shape[0]
        return pcount

    def _identify_fields(self, data_file):
        return self.fields[data_file.filename].keys(), {}

class IOHandlerStreamHexahedral(BaseIOHandler):
    _dataset_type = "stream_hexahedral"
    _vector_fields = ("particle_velocity", "particle_position")

    def __init__(self, ds):
        self.fields = ds.stream_handler.fields
        super(IOHandlerStreamHexahedral, self).__init__(ds)

    def _read_fluid_selection(self, chunks, selector, fields, size):
        chunks = list(chunks)
        assert(len(chunks) == 1)
        chunk = chunks[0]
        rv = {}
        for field in fields:
            ftype, fname = field
            rv[field] = np.empty(size, dtype="float64")
        ngrids = sum(len(chunk.objs) for chunk in chunks)
        mylog.debug("Reading %s cells of %s fields in %s blocks",
                    size, [fn for ft, fn in fields], ngrids)
        for field in fields:
            ind = 0
            ftype, fname = field
            for chunk in chunks:
                for g in chunk.objs:
                    ds = self.fields[g.mesh_id].get(field, None)
                    if ds is None:
                        ds = self.fields[g.mesh_id][fname]
                    ind += g.select(selector, ds, rv[field], ind) # caches
        return rv

class IOHandlerStreamOctree(BaseIOHandler):
    _dataset_type = "stream_octree"
    _vector_fields = ("particle_velocity", "particle_position")

    def __init__(self, ds):
        self.fields = ds.stream_handler.fields
        super(IOHandlerStreamOctree, self).__init__(ds)

    def _read_fluid_selection(self, chunks, selector, fields, size):
        rv = {}
        ind = 0
        chunks = list(chunks)
        assert(len(chunks) == 1)
        for chunk in chunks:
            assert(len(chunk.objs) == 1)
            for subset in chunk.objs:
                field_vals = {}
                for field in fields:
                    field_vals[field] = self.fields[
                        subset.domain_id - subset._domain_offset][field]
                subset.fill(field_vals, rv, selector, ind)
        return rv


class IOHandlerStreamUnstructured(BaseIOHandler):
    _dataset_type = "stream_unstructured"

    def __init__(self, ds):
        self.fields = ds.stream_handler.fields
        super(IOHandlerStreamUnstructured, self).__init__(ds)

    def _read_fluid_selection(self, chunks, selector, fields, size):
        chunks = list(chunks)
        chunk = chunks[0]
        mesh_id = chunk.objs[0].mesh_id
        rv = {}
        for field in fields:
            if field in self.ds._node_fields:
                nodes_per_element = self.fields[mesh_id][field].shape[1]
                rv[field] = np.empty((size, nodes_per_element), dtype="float64")
            else:
                rv[field] = np.empty(size, dtype="float64")
        ngrids = sum(len(chunk.objs) for chunk in chunks)
        mylog.debug("Reading %s cells of %s fields in %s blocks",
                    size, [fname for ftype, fname in fields], ngrids)
        for field in fields:
            ind = 0
            ftype, fname = field
            for chunk in chunks:
                for g in chunk.objs:
                    ds = self.fields[g.mesh_id].get(field, None)
                    if ds is None:
                        ds = self.fields[g.mesh_id][fname]
                    ind += g.select(selector, ds, rv[field], ind) # caches
        return rv

