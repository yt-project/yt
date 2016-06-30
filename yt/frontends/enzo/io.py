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

import random
from contextlib import contextmanager

from yt.utilities.io_handler import \
    BaseIOHandler
from yt.utilities.logger import ytLogger as mylog
from yt.geometry.selection_routines import AlwaysSelector
from yt.extern.six import u, b, iteritems
from yt.utilities.on_demand_imports import _h5py as h5py

import numpy as np


_convert_mass = ("particle_mass","mass")

_particle_position_names = {}

class IOHandlerPackedHDF5(BaseIOHandler):

    _dataset_type = "enzo_packed_3d"
    _base = slice(None)
    _field_dtype = "float64"

    def _read_field_names(self, grid):
        if grid.filename is None: return []
        f = h5py.File(grid.filename, "r")
        try:
            group = f["/Grid%08i" % grid.id]
        except KeyError:
            group = f
        fields = []
        dtypes = set([])
        add_io = "io" in grid.ds.particle_types
        for name, v in iteritems(group):
            # NOTE: This won't work with 1D datasets or references.
            # For all versions of Enzo I know about, we can assume all floats
            # are of the same size.  So, let's grab one.
            if not hasattr(v, "shape") or v.dtype == "O":
                continue
            elif len(v.dims) == 1:
                if grid.ds.dimensionality == 1:
                    fields.append( ("enzo", str(name)) )
                elif add_io:
                    fields.append( ("io", str(name)) )
            else:
                fields.append( ("enzo", str(name)) )
                dtypes.add(v.dtype)

        if len(dtypes) == 1:
            # Now, if everything we saw was the same dtype, we can go ahead and
            # set it here.  We do this because it is a HUGE savings for 32 bit
            # floats, since our numpy copying/casting is way faster than
            # h5py's, for some reason I don't understand.  This does *not* need
            # to be correct -- it will get fixed later -- it just needs to be
            # okay for now.
            self._field_dtype = list(dtypes)[0]
        f.close()
        return fields

    @property
    def _read_exception(self):
        return (KeyError,)

    def _read_particle_coords(self, chunks, ptf):
        chunks = list(chunks)
        for chunk in chunks: # These should be organized by grid filename
            f = None
            for g in chunk.objs:
                if g.filename is None: continue
                if f is None:
                    #print "Opening (count) %s" % g.filename
                    f = h5py.File(g.filename, "r")
                nap = sum(g.NumberOfActiveParticles.values())
                if g.NumberOfParticles == 0 and nap == 0:
                    continue
                ds = f.get("/Grid%08i" % g.id)
                for ptype, field_list in sorted(ptf.items()):
                    if ptype != "io":
                        if g.NumberOfActiveParticles[ptype] == 0: continue
                        pds = ds.get("Particles/%s" % ptype)
                    else:
                        pds = ds
                    pn = _particle_position_names.get(ptype,
                            r"particle_position_%s")
                    x, y, z = (np.asarray(pds.get(pn % ax).value, dtype="=f8")
                               for ax in 'xyz')
                    for field in field_list:
                        if np.asarray(pds[field]).ndim > 1:
                            self._array_fields[field] = pds[field].shape
                    yield ptype, (x, y, z)
            if f: f.close()

    def _read_particle_fields(self, chunks, ptf, selector):
        chunks = list(chunks)
        for chunk in chunks: # These should be organized by grid filename
            f = None
            for g in chunk.objs:
                if g.filename is None: continue
                if f is None:
                    #print "Opening (read) %s" % g.filename
                    f = h5py.File(g.filename, "r")
                nap = sum(g.NumberOfActiveParticles.values())
                if g.NumberOfParticles == 0 and nap == 0:
                    continue
                ds = f.get("/Grid%08i" % g.id)
                for ptype, field_list in sorted(ptf.items()):
                    if ptype != "io":
                        if g.NumberOfActiveParticles[ptype] == 0: continue
                        pds = ds.get("Particles/%s" % ptype)
                    else:
                        pds = ds
                    pn = _particle_position_names.get(ptype,
                            r"particle_position_%s")
                    x, y, z = (np.asarray(pds.get(pn % ax).value, dtype="=f8")
                               for ax in 'xyz')
                    mask = selector.select_points(x, y, z, 0.0)
                    if mask is None: continue
                    for field in field_list:
                        data = np.asarray(pds.get(field).value, "=f8")
                        if field in _convert_mass:
                            data *= g.dds.prod(dtype="f8")
                        yield (ptype, field), data[mask]
            if f: f.close()

    def _read_fluid_selection(self, chunks, selector, fields, size):
        rv = {}
        # Now we have to do something unpleasant
        chunks = list(chunks)
        if selector.__class__.__name__ == "GridSelector":
            if not (len(chunks) == len(chunks[0].objs) == 1):
                raise RuntimeError
            g = chunks[0].objs[0]
            f = h5py.File(u(g.filename), 'r')
            if g.id in self._cached_fields:
                gf = self._cached_fields[g.id]
                rv.update(gf)
            if len(rv) == len(fields): return rv
            gds = f.get("/Grid%08i" % g.id)
            for field in fields:
                if field in rv:
                    self._hits += 1
                    continue
                self._misses += 1
                ftype, fname = field
                if fname in gds:
                    rv[(ftype, fname)] = gds.get(fname).value.swapaxes(0, -1)
                else:
                    rv[(ftype, fname)] = np.zeros(g.ActiveDimensions)
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
        h5_type = self._field_dtype
        for chunk in chunks:
            fid = None
            for g in chunk.objs:
                if g.filename is None: continue
                if fid is None:
                    fid = h5py.h5f.open(b(g.filename), h5py.h5f.ACC_RDONLY)
                gf = self._cached_fields.get(g.id, {})
                data = np.empty(g.ActiveDimensions[::-1], dtype=h5_type)
                data_view = data.swapaxes(0, -1)
                nd = 0
                for field in fields:
                    if field in gf:
                        nd = g.select(selector, gf[field], rv[field], ind)
                        self._hits += 1
                        continue
                    self._misses += 1
                    ftype, fname = field
                    try:
                        node = "/Grid%08i/%s" % (g.id, fname)
                        dg = h5py.h5d.open(fid, b(node))
                    except KeyError:
                        if fname == "Dark_Matter_Density": continue
                        raise
                    dg.read(h5py.h5s.ALL, h5py.h5s.ALL, data)
                    if self._cache_on:
                        self._cached_fields.setdefault(g.id, {})
                        # Copy because it's a view into an empty temp array
                        self._cached_fields[g.id][field] = data_view.copy()
                    nd = g.select(selector, data_view, rv[field], ind) # caches
                ind += nd
            if fid: fid.close()
        return rv

    @contextmanager
    def preload(self, chunk, fields, max_size):
        if len(fields) == 0:
            yield self
            return
        old_cache_on = self._cache_on
        old_cached_fields = self._cached_fields
        self._cached_fields = cf = {}
        self._cache_on = True
        for gid in old_cached_fields:
            # Will not copy numpy arrays, which is good!
            cf[gid] = old_cached_fields[gid].copy() 
        self._hits = self._misses = 0
        self._cached_fields = self._read_chunk_data(chunk, fields)
        mylog.debug("(1st) Hits = % 10i Misses = % 10i",
            self._hits, self._misses)
        self._hits = self._misses = 0
        yield self
        mylog.debug("(2nd) Hits = % 10i Misses = % 10i",
            self._hits, self._misses)
        self._cached_fields = old_cached_fields
        self._cache_on = old_cache_on
        # Randomly remove some grids from the cache.  Note that we're doing
        # this on a grid basis, not a field basis.  Performance will be
        # slightly non-deterministic as a result of this, but it should roughly
        # be statistically alright, assuming (as we do) that this will get
        # called during largely unbalanced stuff.
        if len(self._cached_fields) > max_size:
            to_remove = random.sample(self._cached_fields.keys(),
                len(self._cached_fields) - max_size)
            mylog.debug("Purging from cache %s", len(to_remove))
            for k in to_remove:
                self._cached_fields.pop(k)
        else:
            mylog.warning("Cache size % 10i (max % 10i)",
                len(self._cached_fields), max_size)

    def _read_chunk_data(self, chunk, fields):
        fid = fn = None
        rv = {}
        mylog.debug("Preloading fields %s", fields)
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
        h5_type = self._field_dtype
        for g in chunk.objs:
            rv[g.id] = gf = {}
            if g.id in self._cached_fields:
                rv[g.id].update(self._cached_fields[g.id])
            if g.filename is None: continue
            elif g.filename != fn:
                if fid is not None: fid.close()
                fid = None
            if fid is None:
                fid = h5py.h5f.open(b(g.filename), h5py.h5f.ACC_RDONLY)
                fn = g.filename
            data = np.empty(g.ActiveDimensions[::-1], dtype=h5_type)
            data_view = data.swapaxes(0, -1)
            for field in fluid_fields:
                if field in gf:
                    self._hits += 1
                    continue
                self._misses += 1
                ftype, fname = field
                try:
                    node = "/Grid%08i/%s" % (g.id, fname)
                    dg = h5py.h5d.open(fid, b(node))
                except KeyError:
                    if fname == "Dark_Matter_Density": continue
                    raise
                dg.read(h5py.h5s.ALL, h5py.h5s.ALL, data)
                gf[field] = data_view.copy()
        if fid: fid.close()
        if self._cache_on:
            for gid in rv:
                self._cached_fields.setdefault(gid, {})
                self._cached_fields[gid].update(rv[gid])
        return rv

class IOHandlerPackedHDF5GhostZones(IOHandlerPackedHDF5):
    _dataset_type = "enzo_packed_3d_gz"

    def __init__(self, *args, **kwargs):
        super(IOHandlerPackedHDF5GhostZones, self).__init__(*args, **kwargs)
        NGZ = self.ds.parameters.get("NumberOfGhostZones", 3)
        self._base = (slice(NGZ, -NGZ),
                      slice(NGZ, -NGZ),
                      slice(NGZ, -NGZ))

    def _read_raw_data_set(self, grid, field):
        f = h5py.File(grid.filename, "r")
        ds = f["/Grid%08i/%s" % (grid.id, field)][:].swapaxes(0,2)
        f.close()
        return ds

class IOHandlerInMemory(BaseIOHandler):

    _dataset_type = "enzo_inline"

    def __init__(self, ds, ghost_zones=3):
        self.ds = ds
        import enzo
        self.enzo = enzo
        self.grids_in_memory = enzo.grid_data
        self.old_grids_in_memory = enzo.old_grid_data
        self.my_slice = (slice(ghost_zones,-ghost_zones),
                      slice(ghost_zones,-ghost_zones),
                      slice(ghost_zones,-ghost_zones))
        BaseIOHandler.__init__(self, ds)

    def _read_field_names(self, grid):
        fields = []
        add_io = "io" in grid.ds.particle_types
        for name, v in self.grids_in_memory[grid.id].items():
            # NOTE: This won't work with 1D datasets or references.
            if not hasattr(v, "shape") or v.dtype == "O":
                continue
            elif v.ndim == 1:
                if grid.ds.dimensionality == 1:
                    fields.append( ("enzo", str(name)) )
                elif add_io:
                    fields.append( ("io", str(name)) )
            else:
                fields.append( ("enzo", str(name)) )
        return fields

    def _read_fluid_selection(self, chunks, selector, fields, size):
        rv = {}
        # Now we have to do something unpleasant
        chunks = list(chunks)
        if selector.__class__.__name__ == "GridSelector":
            if not (len(chunks) == len(chunks[0].objs) == 1):
                raise RuntimeError
            g = chunks[0].objs[0]
            for ftype, fname in fields:
                rv[(ftype, fname)] = \
                    self.grids_in_memory[g.id][fname].swapaxes(0, 2)
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
            for g in chunk.objs:
                # We want a *hard error* here.
                #if g.id not in self.grids_in_memory: continue
                for field in fields:
                    ftype, fname = field
                    data_view = self.grids_in_memory[g.id][fname][self.my_slice].swapaxes(0,2)
                    nd = g.select(selector, data_view, rv[field], ind)
                ind += nd
        assert(ind == fsize)
        return rv

    def _read_particle_coords(self, chunks, ptf):
        chunks = list(chunks)
        for chunk in chunks: # These should be organized by grid filename
            for g in chunk.objs:
                if g.id not in self.grids_in_memory: continue
                nap = sum(g.NumberOfActiveParticles.values())
                if g.NumberOfParticles == 0 and nap == 0: continue
                for ptype, field_list in sorted(ptf.items()):
                    x, y, z = self.grids_in_memory[g.id]['particle_position_x'], \
                                        self.grids_in_memory[g.id]['particle_position_y'], \
                                        self.grids_in_memory[g.id]['particle_position_z']
                    yield ptype, (x, y, z)

    def _read_particle_fields(self, chunks, ptf, selector):
        chunks = list(chunks)
        for chunk in chunks: # These should be organized by grid filename
            for g in chunk.objs:
                if g.id not in self.grids_in_memory: continue
                nap = sum(g.NumberOfActiveParticles.values())
                if g.NumberOfParticles == 0 and nap == 0: continue
                for ptype, field_list in sorted(ptf.items()):
                    x, y, z = self.grids_in_memory[g.id]['particle_position_x'], \
                                        self.grids_in_memory[g.id]['particle_position_y'], \
                                        self.grids_in_memory[g.id]['particle_position_z']
                    mask = selector.select_points(x, y, z, 0.0)
                    if mask is None: continue
                    for field in field_list:
                        data = self.grids_in_memory[g.id][field]
                        if field in _convert_mass:
                            data = data * g.dds.prod(dtype="f8")
                        yield (ptype, field), data[mask]

class IOHandlerPacked2D(IOHandlerPackedHDF5):

    _dataset_type = "enzo_packed_2d"
    _particle_reader = False

    def _read_data_set(self, grid, field):
        f = h5py.File(grid.filename, "r")
        ds = f["/Grid%08i/%s" % (grid.id, field)][:]
        f.close()
        return ds.transpose()[:,:,None]

    def modify(self, field):
        pass

    def _read_fluid_selection(self, chunks, selector, fields, size):
        rv = {}
        # Now we have to do something unpleasant
        chunks = list(chunks)
        if selector.__class__.__name__ == "GridSelector":
            if not (len(chunks) == len(chunks[0].objs) == 1):
                raise RuntimeError
            g = chunks[0].objs[0]
            f = h5py.File(g.filename, 'r')
            gds = f.get("/Grid%08i" % g.id)
            for ftype, fname in fields:
                rv[(ftype, fname)] = np.atleast_3d(
                    gds.get(fname).value.transpose())
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
                if f is None:
                    #print "Opening (count) %s" % g.filename
                    f = h5py.File(g.filename, "r")
                gds = f.get("/Grid%08i" % g.id)
                if gds is None:
                    gds = f
                for field in fields:
                    ftype, fname = field
                    ds = np.atleast_3d(gds.get(fname).value.transpose())
                    nd = g.select(selector, ds, rv[field], ind) # caches
                ind += nd
            f.close()
        return rv

class IOHandlerPacked1D(IOHandlerPackedHDF5):

    _dataset_type = "enzo_packed_1d"
    _particle_reader = False

    def _read_data_set(self, grid, field):
        f = h5py.File(grid.filename, "r")
        ds = f["/Grid%08i/%s" % (grid.id, field)][:]
        f.close()
        return ds.transpose()[:,None,None]

    def modify(self, field):
        pass
