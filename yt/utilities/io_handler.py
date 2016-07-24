"""
The data-file handling functions



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from collections import defaultdict
from contextlib import contextmanager

import os
from yt.utilities.on_demand_imports import _h5py as h5py
import numpy as np
from yt.extern.six import add_metaclass

_axis_ids = {0:2,1:1,2:0}

io_registry = {}

class RegisteredIOHandler(type):
    def __init__(cls, name, b, d):
        type.__init__(cls, name, b, d)
        if hasattr(cls, "_dataset_type"):
            io_registry[cls._dataset_type] = cls

@add_metaclass(RegisteredIOHandler)
class BaseIOHandler(object):
    _vector_fields = ()
    _dataset_type = None
    _particle_reader = False
    _cache_on = False
    _misses = 0
    _hits = 0

    def __init__(self, ds):
        self.queue = defaultdict(dict)
        self.ds = ds
        self._last_selector_id = None
        self._last_selector_counts = None
        self._array_fields = {}
        self._cached_fields = {}
        # Make sure _vector_fields is a dict of fields and their dimension
        # and assume all non-specified vector fields are 3D
        if not isinstance(self._vector_fields, dict):
            self._vector_fields = dict((field, 3) for field in self._vector_fields)

    # We need a function for reading a list of sets
    # and a function for *popping* from a queue all the appropriate sets

    @contextmanager
    def preload(self, chunk, fields, max_size):
        yield self

    def pop(self, grid, field):
        if grid.id in self.queue and field in self.queue[grid.id]:
            return self.modify(self.queue[grid.id].pop(field))
        else:
            # We only read the one set and do not store it if it isn't pre-loaded
            return self._read_data_set(grid, field)

    def peek(self, grid, field):
        return self.queue[grid.id].get(field, None)

    def push(self, grid, field, data):
        if grid.id in self.queue and field in self.queue[grid.id]:
            raise ValueError
        self.queue[grid][field] = data

    def _field_in_backup(self, grid, backup_file, field_name):
        if os.path.exists(backup_file):
            fhandle = h5py.File(backup_file, 'r')
            g = fhandle["data"]
            grid_group = g["grid_%010i" % (grid.id - grid._id_offset)]
            if field_name in grid_group:
                return_val = True
            else:
                return_val = False
            fhandle.close()
            return return_val
        else:
            return False
            
    def _read_data_set(self, grid, field):
        # check backup file first. if field not found,
        # call frontend-specific io method
        backup_filename = grid.ds.backup_filename
        if not grid.ds.read_from_backup:
            return self._read_data(grid, field)
        elif self._field_in_backup(grid, backup_filename, field):
            fhandle = h5py.File(backup_filename, 'r')
            g = fhandle["data"]
            grid_group = g["grid_%010i" % (grid.id - grid._id_offset)]
            data = grid_group[field][:]
            fhandle.close()
            return data
        else:
            return self._read_data(grid, field)
                
    # Now we define our interface
    def _read_data(self, grid, field):
        pass

    def _read_data_slice(self, grid, field, axis, coord):
        sl = [slice(None), slice(None), slice(None)]
        sl[axis] = slice(coord, coord + 1)
        tr = self._read_data_set(grid, field)[sl]
        if tr.dtype == "float32": tr = tr.astype("float64")
        return tr

    def _read_field_names(self, grid):
        pass

    @property
    def _read_exception(self):
        return None

    def _read_chunk_data(self, chunk, fields):
        return {}

    def _count_particles_chunks(self, chunks, ptf, selector):
        psize = defaultdict(lambda: 0) # COUNT PTYPES ON DISK
        for ptype, (x, y, z) in self._read_particle_coords(chunks, ptf):
            psize[ptype] += selector.count_points(x, y, z, 0.0)
        return dict(psize.items())

    def _read_particle_selection(self, chunks, selector, fields):
        rv = {}
        ind = {}
        # We first need a set of masks for each particle type
        ptf = defaultdict(list)        # ON-DISK TO READ
        fsize = defaultdict(lambda: 0) # COUNT RV
        field_maps = defaultdict(list) # ptypes -> fields
        chunks = list(chunks)
        unions = self.ds.particle_unions
        # What we need is a mapping from particle types to return types
        for field in fields:
            ftype, fname = field
            fsize[field] = 0
            # We should add a check for p.fparticle_unions or something here
            if ftype in unions:
                for pt in unions[ftype]:
                    ptf[pt].append(fname)
                    field_maps[pt, fname].append(field)
            else:
                ptf[ftype].append(fname)
                field_maps[field].append(field)
        # We can't hash chunks, but otherwise this is a neat idea.
        # Now we have our full listing.
        # Here, ptype_map means which particles contribute to a given type.
        # And ptf is the actual fields from disk to read.
        psize = self._count_particles_chunks(chunks, ptf, selector)
        # Now we allocate
        # ptf, remember, is our mapping of what we want to read
        #for ptype in ptf:
        for field in fields:
            if field[0] in unions:
                for pt in unions[field[0]]:
                    fsize[field] += psize.get(pt, 0)
            else:
                fsize[field] += psize.get(field[0], 0)
        for field in fields:
            if field[1] in self._vector_fields:
                shape = (fsize[field], self._vector_fields[field[1]])
            elif field[1] in self._array_fields:
                shape = (fsize[field],)+self._array_fields[field[1]]
            else:
                shape = (fsize[field], )
            rv[field] = np.empty(shape, dtype="float64")
            ind[field] = 0
        # Now we read.
        for field_r, vals in self._read_particle_fields(chunks, ptf, selector):
            # Note that we now need to check the mappings
            for field_f in field_maps[field_r]:
                my_ind = ind[field_f]
                #mylog.debug("Filling %s from %s to %s with %s",
                #    field_f, my_ind, my_ind+vals.shape[0], field_r)
                rv[field_f][my_ind:my_ind + vals.shape[0],...] = vals
                ind[field_f] += vals.shape[0]
        # Now we need to truncate all our fields, since we allow for
        # over-estimating.
        for field_f in ind:
            rv[field_f] = rv[field_f][:ind[field_f]]
        return rv

class IOHandlerExtracted(BaseIOHandler):

    _dataset_type = 'extracted'

    def _read_data_set(self, grid, field):
        return (grid.base_grid[field] / grid.base_grid.convert(field))

    def _read_data_slice(self, grid, field, axis, coord):
        sl = [slice(None), slice(None), slice(None)]
        sl[axis] = slice(coord, coord + 1)
        return grid.base_grid[field][tuple(sl)] / grid.base_grid.convert(field)
