"""
YTData data-file handling function




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2015, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from yt.extern.six import \
    u
from yt.funcs import \
    mylog
from yt.utilities.exceptions import \
    YTDomainOverflow
from yt.utilities.io_handler import \
    BaseIOHandler
from yt.utilities.lib.geometry_utils import \
    compute_morton
from yt.utilities.on_demand_imports import \
    _h5py as h5py

class IOHandlerYTNonspatialhdf5(BaseIOHandler):
    _dataset_type = "ytnonspatialhdf5"
    _base = slice(None)
    _field_dtype = "float64"

    def _read_fluid_selection(self, g, selector, fields):
        rv = {}
        if selector.__class__.__name__ == "GridSelector":
            if g.id in self._cached_fields:
                gf = self._cached_fields[g.id]
                rv.update(gf)
            if len(rv) == len(fields): return rv
            f = h5py.File(u(g.filename), "r")
            gds = f["data"]
            for field in fields:
                if field in rv:
                    self._hits += 1
                    continue
                self._misses += 1
                ftype, fname = field
                rv[(ftype, fname)] = gds[fname].value
            if self._cache_on:
                for gid in rv:
                    self._cached_fields.setdefault(gid, {})
                    self._cached_fields[gid].update(rv[gid])
            f.close()
            return rv
        else:
            raise RuntimeError(
                "Geometric selection not supported for non-spatial datasets.")

class IOHandlerYTGridHDF5(BaseIOHandler):
    _dataset_type = "ytgridhdf5"
    _base = slice(None)
    _field_dtype = "float64"

    def _read_fluid_selection(self, chunks, selector, fields, size):
        rv = {}
        # Now we have to do something unpleasant
        chunks = list(chunks)
        if selector.__class__.__name__ == "GridSelector":
            if not (len(chunks) == len(chunks[0].objs) == 1):
                raise RuntimeError
            g = chunks[0].objs[0]
            if g.id in self._cached_fields:
                gf = self._cached_fields[g.id]
                rv.update(gf)
            if len(rv) == len(fields): return rv
            f = h5py.File(u(g.filename), "r")
            gds = f[self.ds.default_fluid_type]
            for field in fields:
                if field in rv:
                    self._hits += 1
                    continue
                self._misses += 1
                ftype, fname = field
                rv[(ftype, fname)] = gds[fname].value
            if self._cache_on:
                for gid in rv:
                    self._cached_fields.setdefault(gid, {})
                    self._cached_fields[gid].update(rv[gid])
            f.close()
            return rv
        if size is None:
            size = sum((g.count(selector) for chunk in chunks
                        for g in chunk.objs))
        for field in fields:
            ftype, fname = field
            fsize = size
            rv[field] = np.empty(fsize, dtype="float64")
        ng = sum(len(c.objs) for c in chunks)
        mylog.debug("Reading %s cells of %s fields in %s grids",
                   size, [f2 for f1, f2 in fields], ng)
        ind = 0
        for chunk in chunks:
            f = None
            for g in chunk.objs:
                if g.filename is None: continue
                if f is None:
                    f = h5py.File(g.filename, "r")
                gf = self._cached_fields.get(g.id, {})
                nd = 0
                for field in fields:
                    if field in gf:
                        nd = g.select(selector, gf[field], rv[field], ind)
                        self._hits += 1
                        continue
                    self._misses += 1
                    ftype, fname = field
                    # add extra dimensions to make data 3D
                    data = f[ftype][fname].value.astype(self._field_dtype)
                    for dim in range(len(data.shape), 3):
                        data = np.expand_dims(data, dim)
                    if self._cache_on:
                        self._cached_fields.setdefault(g.id, {})
                        self._cached_fields[g.id][field] = data
                    nd = g.select(selector, data, rv[field], ind) # caches
                ind += nd
            if f: f.close()
        return rv

    def _read_particle_coords(self, chunks, ptf):
        chunks = list(chunks)
        for chunk in chunks:
            f = None
            for g in chunk.objs:
                if g.filename is None: continue
                if f is None:
                    f = h5py.File(g.filename, "r")
                if g.NumberOfParticles == 0:
                    continue
                for ptype, field_list in sorted(ptf.items()):
                    pn = "particle_position_%s"
                    x, y, z = (np.asarray(f[ptype][pn % ax].value, dtype="=f8")
                               for ax in 'xyz')
                    for field in field_list:
                        if np.asarray(f[ptype][field]).ndim > 1:
                            self._array_fields[field] = f[ptype][field].shape
                    yield ptype, (x, y, z)
            if f: f.close()

    def _read_particle_fields(self, chunks, ptf, selector):
        chunks = list(chunks)
        for chunk in chunks: # These should be organized by grid filename
            f = None
            for g in chunk.objs:
                if g.filename is None: continue
                if f is None:
                    f = h5py.File(g.filename, "r")
                if g.NumberOfParticles == 0:
                    continue
                for ptype, field_list in sorted(ptf.items()):
                    pn = "particle_position_%s"
                    x, y, z = (np.asarray(f[ptype][pn % ax].value, dtype="=f8")
                               for ax in 'xyz')
                    mask = selector.select_points(x, y, z, 0.0)
                    if mask is None: continue
                    for field in field_list:
                        data = np.asarray(f[ptype][field].value, "=f8")
                        yield (ptype, field), data[mask]
            if f: f.close()

class IOHandlerYTDataContainerHDF5(BaseIOHandler):
    _dataset_type = "ytdatacontainer_hdf5"

    def _read_fluid_selection(self, chunks, selector, fields, size):
        raise NotImplementedError

    def _read_particle_coords(self, chunks, ptf):
        # This will read chunks and yield the results.
        chunks = list(chunks)
        data_files = set([])
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        for data_file in sorted(data_files):
            with h5py.File(data_file.filename, "r") as f:
                for ptype, field_list in sorted(ptf.items()):
                    pcount = data_file.total_particles[ptype]
                    if pcount == 0: continue
                    x = _get_position_array(ptype, f, "x")
                    y = _get_position_array(ptype, f, "y")
                    z = _get_position_array(ptype, f, "z")
                    yield ptype, (x, y, z)

    def _read_particle_fields(self, chunks, ptf, selector):
        # Now we have all the sizes, and we can allocate
        chunks = list(chunks)
        data_files = set([])
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        for data_file in sorted(data_files):
            with h5py.File(data_file.filename, "r") as f:
                for ptype, field_list in sorted(ptf.items()):
                    x = _get_position_array(ptype, f, "x")
                    y = _get_position_array(ptype, f, "y")
                    z = _get_position_array(ptype, f, "z")
                    mask = selector.select_points(x, y, z, 0.0)
                    del x, y, z
                    if mask is None: continue
                    for field in field_list:
                        data = f[ptype][field][mask].astype("float64")
                        yield (ptype, field), data

    def _initialize_index(self, data_file, regions):
        all_count = self._count_particles(data_file)
        pcount = sum(all_count.values())
        morton = np.empty(pcount, dtype='uint64')
        mylog.debug("Initializing index % 5i (% 7i particles)",
                    data_file.file_id, pcount)
        ind = 0
        with h5py.File(data_file.filename, "r") as f:
            for ptype in all_count:
                if ptype not in f or all_count[ptype] == 0: continue
                pos = np.empty((all_count[ptype], 3), dtype="float64")
                pos = data_file.ds.arr(pos, "code_length")
                if ptype == "grid":
                    dx = f["grid"]["dx"].value.min()
                else:
                    dx = 2. * np.finfo(f[ptype]["particle_position_x"].dtype).eps
                dx = self.ds.quan(dx, "code_length")
                pos[:,0] = _get_position_array(ptype, f, "x")
                pos[:,1] = _get_position_array(ptype, f, "y")
                pos[:,2] = _get_position_array(ptype, f, "z")
                # These are 32 bit numbers, so we give a little lee-way.
                # Otherwise, for big sets of particles, we often will bump into the
                # domain edges.  This helps alleviate that.
                np.clip(pos, self.ds.domain_left_edge + dx,
                             self.ds.domain_right_edge - dx, pos)
                if np.any(pos.min(axis=0) < self.ds.domain_left_edge) or \
                   np.any(pos.max(axis=0) > self.ds.domain_right_edge):
                    raise YTDomainOverflow(pos.min(axis=0),
                                           pos.max(axis=0),
                                           self.ds.domain_left_edge,
                                           self.ds.domain_right_edge)
                regions.add_data_file(pos, data_file.file_id)
                morton[ind:ind+pos.shape[0]] = compute_morton(
                    pos[:,0], pos[:,1], pos[:,2],
                    data_file.ds.domain_left_edge,
                    data_file.ds.domain_right_edge)
                ind += pos.shape[0]
        return morton

    def _count_particles(self, data_file):
        return self.ds.num_particles

    def _identify_fields(self, data_file):
        fields = []
        units = {}
        with h5py.File(data_file.filename, "r") as f:
            for ptype in f:
                fields.extend([(ptype, str(field)) for field in f[ptype]])
                units.update(dict([((ptype, str(field)), 
                                    f[ptype][field].attrs["units"])
                                   for field in f[ptype]]))
        return fields, units

class IOHandlerYTSpatialPlotHDF5(IOHandlerYTDataContainerHDF5):
    _dataset_type = "ytspatialplot_hdf5"

    def _read_particle_coords(self, chunks, ptf):
        # This will read chunks and yield the results.
        chunks = list(chunks)
        data_files = set([])
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        for data_file in sorted(data_files):
            with h5py.File(data_file.filename, "r") as f:
                for ptype, field_list in sorted(ptf.items()):
                    pcount = data_file.total_particles[ptype]
                    if pcount == 0: continue
                    x = _get_position_array(ptype, f, "px")
                    y = _get_position_array(ptype, f, "py")
                    z = np.zeros(x.size, dtype="float64") + \
                      self.ds.domain_left_edge[2].in_cgs().d
                    yield ptype, (x, y, z)

    def _read_particle_fields(self, chunks, ptf, selector):
        # Now we have all the sizes, and we can allocate
        chunks = list(chunks)
        data_files = set([])
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        for data_file in sorted(data_files):
            all_count = self._count_particles(data_file)
            with h5py.File(data_file.filename, "r") as f:
                for ptype, field_list in sorted(ptf.items()):
                    x = _get_position_array(ptype, f, "px")
                    y = _get_position_array(ptype, f, "py")
                    z = np.zeros(all_count[ptype], dtype="float64") + \
                      self.ds.domain_left_edge[2].in_cgs().d
                    mask = selector.select_points(x, y, z, 0.0)
                    del x, y, z
                    if mask is None: continue
                    for field in field_list:
                        data = f[ptype][field][mask].astype("float64")
                        yield (ptype, field), data

    def _initialize_index(self, data_file, regions):
        all_count = self._count_particles(data_file)
        pcount = sum(all_count.values())
        morton = np.empty(pcount, dtype='uint64')
        mylog.debug("Initializing index % 5i (% 7i particles)",
                    data_file.file_id, pcount)
        ind = 0
        with h5py.File(data_file.filename, "r") as f:
            for ptype in all_count:
                if ptype not in f or all_count[ptype] == 0: continue
                pos = np.empty((all_count[ptype], 3), dtype="float64")
                pos = data_file.ds.arr(pos, "code_length")
                if ptype == "grid":
                    dx = f["grid"]["pdx"].value.min()
                else:
                    raise NotImplementedError
                dx = self.ds.quan(dx, "code_length")
                pos[:,0] = _get_position_array(ptype, f, "px")
                pos[:,1] = _get_position_array(ptype, f, "py")
                pos[:,2] = np.zeros(all_count[ptype], dtype="float64") + \
                  self.ds.domain_left_edge[2].in_cgs().d
                # These are 32 bit numbers, so we give a little lee-way.
                # Otherwise, for big sets of particles, we often will bump into the
                # domain edges.  This helps alleviate that.
                np.clip(pos, self.ds.domain_left_edge + dx,
                             self.ds.domain_right_edge - dx, pos)
                if np.any(pos.min(axis=0) < self.ds.domain_left_edge) or \
                   np.any(pos.max(axis=0) > self.ds.domain_right_edge):
                    raise YTDomainOverflow(pos.min(axis=0),
                                           pos.max(axis=0),
                                           self.ds.domain_left_edge,
                                           self.ds.domain_right_edge)
                regions.add_data_file(pos, data_file.file_id)
                morton[ind:ind+pos.shape[0]] = compute_morton(
                    pos[:,0], pos[:,1], pos[:,2],
                    data_file.ds.domain_left_edge,
                    data_file.ds.domain_right_edge)
                ind += pos.shape[0]
        return morton

def _get_position_array(ptype, f, ax):
    if ptype == "grid":
        pos_name = ""
    else:
        pos_name = "particle_position_"
    return f[ptype][pos_name + ax].value.astype("float64")
