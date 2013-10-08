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

import exceptions
import os

from yt.utilities.io_handler import \
    BaseIOHandler, _axis_ids
from yt.utilities.logger import ytLogger as mylog
from yt.geometry.selection_routines import mask_fill
import h5py

import numpy as np
from yt.funcs import *

_convert_mass = ("particle_mass","mass")

_particle_position_names = {}

class IOHandlerPackedHDF5(BaseIOHandler):

    _data_style = "enzo_packed_3d"
    _base = slice(None)

    def _read_field_names(self, grid):
        f = h5py.File(grid.filename, "r")
        group = f["/Grid%08i" % grid.id]
        fields = []
        for name, v in group.iteritems():
            # NOTE: This won't work with 1D datasets.
            if not hasattr(v, "shape"):
                continue
            elif len(v.dims) == 1:
                fields.append( ("io", str(name)) )
            else:
                fields.append( ("gas", str(name)) )
        f.close()
        return fields

    @property
    def _read_exception(self):
        return (exceptions.KeyError,)

    def _read_particle_coords(self, chunks, ptf):
        chunks = list(chunks)
        for chunk in chunks: # These should be organized by grid filename
            f = None
            for g in chunk.objs:
                if f is None:
                    #print "Opening (count) %s" % g.filename
                    f = h5py.File(g.filename, "r")
                if g.NumberOfParticles == 0: continue
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
                    yield ptype, (x, y, z)
            f.close()

    def _read_particle_fields(self, chunks, ptf, selector):
        chunks = list(chunks)
        for chunk in chunks: # These should be organized by grid filename
            f = None
            for g in chunk.objs:
                if f is None:
                    #print "Opening (read) %s" % g.filename
                    f = h5py.File(g.filename, "r")
                if g.NumberOfParticles == 0: continue
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
                    mask = selector.select_points(x, y, z)
                    if mask is None: continue
                    for field in field_list:
                        data = np.asarray(pds.get(field).value, "=f8")
                        if field in _convert_mass:
                            data *= g.dds.prod(dtype="f8")
                        yield (ptype, field), data
            f.close()

    def _read_grid_chunk(self, chunks, fields):
        sets = [fname for ftype, fname in fields]
        g = chunks[0].objs[0]
        if g.filename is None:
            return {}
        rv = hdf5_light_reader.ReadMultipleGrids(
            g.filename, [g.id], sets, "")[g.id]
        for ftype, fname in fields:
            rv[(ftype, fname)] = rv.pop(fname).swapaxes(0,2)
        return rv

    def _read_chunk_data(self, chunk, fields, suffix = ""):
        data = {}
        grids_by_file = defaultdict(list)
        for g in chunk.objs:
            if g.filename is None:
                continue
            grids_by_file[g.filename].append(g.id)
        #if len(chunk.objs) == 1 and len(grids_by_file) > 0:
        #    raise RuntimeError
        sets = [fname for ftype, fname in fields]
        for filename in grids_by_file:
            nodes = grids_by_file[filename]
            nodes.sort()
            data.update(hdf5_light_reader.ReadMultipleGrids(
                filename, nodes, sets, suffix))
        return data

    def _read_fluid_selection(self, chunks, selector, fields, size):
        rv = {}
        # Now we have to do something unpleasant
        chunks = list(chunks)
        if selector.__class__.__name__ == "GridSelector":
            if not (len(chunks) == len(chunks[0].objs) == 1):
                raise RuntimeError
            g = chunks[0].objs[0]
            f = h5py.File(g.filename)
            gds = f.get("/Grid%08i" % g.id)
            for ftype, fname in fields:
                rv[(ftype, fname)] = gds.get(fname).value.swapaxes(0,2)
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
            fid = None
            for g in chunk.objs:
                if fid is None:
                    fid = h5py.h5f.open(g.filename, h5py.h5f.ACC_RDONLY)
                data = np.empty(g.ActiveDimensions[::-1], dtype="float64")
                data_view = data.swapaxes(0,2)
                for field in fields:
                    ftype, fname = field
                    dg = h5py.h5d.open(fid, "/Grid%08i/Density" % g.id)
                    dg.read(h5py.h5s.ALL, h5py.h5s.ALL, data)
                    nd = g.select(selector, data_view, rv[field], ind) # caches
                ind += nd
        return rv


class IOHandlerPackedHDF5GhostZones(IOHandlerPackedHDF5):
    _data_style = "enzo_packed_3d_gz"

    def __init__(self, *args, **kwargs):
        super(IOHandlerPackgedHDF5GhostZones, self).__init__(*args, **kwargs)
        NGZ = self.pf.parameters.get("NumberOfGhostZones", 3)
        self._base = (slice(NGZ, -NGZ),
                      slice(NGZ, -NGZ),
                      slice(NGZ, -NGZ))

    def _read_raw_data_set(self, grid, field):
        raise NotImplementedError
        return hdf5_light_reader.ReadData(grid.filename,
                "/Grid%08i/%s" % (grid.id, field))

class IOHandlerInMemory(BaseIOHandler):

    _data_style = "enzo_inline"

    def __init__(self, pf, ghost_zones=3):
        self.pf = pf
        import enzo
        self.enzo = enzo
        self.grids_in_memory = enzo.grid_data
        self.old_grids_in_memory = enzo.old_grid_data
        self.my_slice = (slice(ghost_zones,-ghost_zones),
                      slice(ghost_zones,-ghost_zones),
                      slice(ghost_zones,-ghost_zones))
        BaseIOHandler.__init__(self, pf)

    def _read_data_set(self, grid, field):
        if grid.id not in self.grids_in_memory:
            mylog.error("Was asked for %s but I have %s", grid.id, self.grids_in_memory.keys())
            raise KeyError
        tr = self.grids_in_memory[grid.id][field]
        # If it's particles, we copy.
        if len(tr.shape) == 1: return tr.copy()
        # New in-place unit conversion breaks if we don't copy first
        return tr.swapaxes(0,2)[self.my_slice].copy()
        # We don't do this, because we currently do not interpolate
        coef1 = max((grid.Time - t1)/(grid.Time - t2), 0.0)
        coef2 = 1.0 - coef1
        t1 = enzo.yt_parameter_file["InitialTime"]
        t2 = enzo.hierarchy_information["GridOldTimes"][grid.id]
        return (coef1*self.grids_in_memory[grid.id][field] + \
                coef2*self.old_grids_in_memory[grid.id][field])\
                [self.my_slice]

    def modify(self, field):
        return field.swapaxes(0,2)

    def _read_field_names(self, grid):
        return self.grids_in_memory[grid.id].keys()

    def _read_data_slice(self, grid, field, axis, coord):
        sl = [slice(3,-3), slice(3,-3), slice(3,-3)]
        sl[axis] = slice(coord + 3, coord + 4)
        sl = tuple(reversed(sl))
        tr = self.grids_in_memory[grid.id][field][sl].swapaxes(0,2)
        # In-place unit conversion requires we return a copy
        return tr.copy()

    @property
    def _read_exception(self):
        return KeyError

class IOHandlerPacked2D(IOHandlerPackedHDF5):

    _data_style = "enzo_packed_2d"
    _particle_reader = False

    def _read_data_set(self, grid, field):
        raise NotImplementedError
        return hdf5_light_reader.ReadData(grid.filename,
            "/Grid%08i/%s" % (grid.id, field)).transpose()[:,:,None]

    def modify(self, field):
        pass

    def _read_data_slice(self, grid, field, axis, coord):
        raise NotImplementedError
        t = hdf5_light_reader.ReadData(grid.filename, "/Grid%08i/%s" %
                        (grid.id, field)).transpose()
        return t

    def _read_fluid_selection(self, chunks, selector, fields, size):
        rv = {}
        # Now we have to do something unpleasant
        chunks = list(chunks)
        if selector.__class__.__name__ == "GridSelector":
            if not (len(chunks) == len(chunks[0].objs) == 1):
                raise RuntimeError
            g = chunks[0].objs[0]
            f = h5py.File(g.filename)
            gds = f.get("/Grid%08i" % g.id)
            for ftype, fname in fields:
                rv[(ftype, fname)] = np.atleast_3d(gds.get(fname).value)
            f.close()
            return rv
        if size is None:
            size = sum((g.count(selector) for chunk in chunks
                        for g in chunk.objs))
        if any((ftype != "gas" for ftype, fname in fields)):
            raise NotImplementedError
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
                for field in fields:
                    ftype, fname = field
                    ds = np.atleast_3d(gds.get(fname).value)
                    nd = g.select(selector, ds, rv[field], ind) # caches
                ind += nd
            f.close()
        return rv

class IOHandlerPacked1D(IOHandlerPackedHDF5):

    _data_style = "enzo_packed_1d"
    _particle_reader = False

    def _read_data_set(self, grid, field):
        raise NotImplementedError
        return hdf5_light_reader.ReadData(grid.filename,
            "/Grid%08i/%s" % (grid.id, field)).transpose()[:,None,None]

    def modify(self, field):
        pass

    def _read_data_slice(self, grid, field, axis, coord):
        raise NotImplementedError
        t = hdf5_light_reader.ReadData(grid.filename, "/Grid%08i/%s" %
                        (grid.id, field))
        return t

