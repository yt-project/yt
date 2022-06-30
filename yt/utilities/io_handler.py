import os
import sys
from collections import defaultdict
from contextlib import contextmanager
from functools import _make_key, lru_cache
from typing import DefaultDict, Dict, List, Tuple

if sys.version_info >= (3, 9):
    from collections.abc import Iterator, Mapping
else:
    from typing import Iterator, Mapping

import numpy as np

from yt._typing import ParticleCoordinateTuple
from yt.geometry.selection_routines import GridSelector
from yt.utilities.on_demand_imports import _h5py as h5py
from yt.utilities.parallel_tools.dask_helper import _get_dask_compute, _get_dask_delayed

io_registry = {}

use_caching = 0


def _make_io_key(args, *_args, **kwargs):
    self, obj, field, ctx = args
    # Ignore self because we have a self-specific cache
    return _make_key((obj.id, field), *_args, **kwargs)


class BaseIOHandler:
    _vector_fields: Dict[str, int] = {}
    _dataset_type: str
    _particle_reader = False
    _cache_on = False
    _misses = 0
    _hits = 0
    _dask_enabled = False  # global dask_enabled config option may override

    def __init_subclass__(cls, *args, **kwargs):
        super().__init_subclass__(*args, **kwargs)
        if hasattr(cls, "_dataset_type"):
            io_registry[cls._dataset_type] = cls
        if use_caching and hasattr(cls, "_read_obj_field"):
            cls._read_obj_field = lru_cache(
                maxsize=use_caching, typed=True, make_key=_make_io_key
            )(cls._read_obj_field)

    def __init__(self, ds):
        self.queue = defaultdict(dict)
        self.ds = ds
        self._last_selector_id = None
        self._last_selector_counts = None
        self._array_fields = {}
        self._cached_fields = {}

    # We need a function for reading a list of sets
    # and a function for *popping* from a queue all the appropriate sets
    @contextmanager
    def preload(self, chunk, fields: List[Tuple[str, str]], max_size):
        yield self

    def peek(self, grid, field):
        return self.queue[grid.id].get(field, None)

    def push(self, grid, field, data):
        if grid.id in self.queue and field in self.queue[grid.id]:
            raise ValueError
        self.queue[grid][field] = data

    def _field_in_backup(self, grid, backup_file, field_name):
        if os.path.exists(backup_file):
            fhandle = h5py.File(backup_file, mode="r")
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
        if not os.path.exists(backup_filename):
            return self._read_data(grid, field)
        elif self._field_in_backup(grid, backup_filename, field):
            fhandle = h5py.File(backup_filename, mode="r")
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

    def _read_fluid_selection(
        self, chunks, selector, fields: List[Tuple[str, str]], size
    ) -> Mapping[Tuple[str, str], np.ndarray]:
        # This function has an interesting history.  It previously was mandate
        # to be defined by all of the subclasses.  But, to avoid having to
        # rewrite a whole bunch of IO handlers all at once, and to allow a
        # better abstraction for grid-based frontends, we're now defining it in
        # the base class.
        rv = {}
        nodal_fields = []
        for field in fields:
            finfo = self.ds.field_info[field]
            nodal_flag = finfo.nodal_flag
            if np.any(nodal_flag):
                num_nodes = 2 ** sum(nodal_flag)
                rv[field] = np.empty((size, num_nodes), dtype="=f8")
                nodal_fields.append(field)
            else:
                rv[field] = np.empty(size, dtype="=f8")
        ind = {field: 0 for field in fields}
        for field, obj, data in self.io_iter(chunks, fields):
            if data is None:
                continue
            if isinstance(selector, GridSelector) and field not in nodal_fields:
                ind[field] += data.size
                rv[field] = data.copy()
            else:
                ind[field] += obj.select(selector, data, rv[field], ind[field])
        return rv

    def io_iter(self, chunks, fields: List[Tuple[str, str]]):
        raise NotImplementedError(
            "subclassing Dataset.io_iter this is required in order to use the default "
            "implementation of Dataset._read_fluid_selection. "
            "Custom implementations of the latter may not rely on this method."
        )

    def _read_data_slice(self, grid, field, axis, coord):
        sl = [slice(None), slice(None), slice(None)]
        sl[axis] = slice(coord, coord + 1)
        tr = self._read_data_set(grid, field)[tuple(sl)]
        if tr.dtype == "float32":
            tr = tr.astype("float64")
        return tr

    def _read_field_names(self, grid):
        pass

    @property
    def _read_exception(self):
        return None

    def _read_chunk_data(self, chunk, fields):
        return {}

    def _read_particle_coords(
        self, chunks, ptf: DefaultDict[str, List[str]]
    ) -> Iterator[ParticleCoordinateTuple]:
        # An iterator that yields particle coordinates for each chunk by particle
        # type. Must be implemented by each frontend. Must yield a tuple of
        # (particle type, xyz, hsml) by chunk. If the frontend does not have
        # a smoothing length, yield (particle type, xyz, 0.0)
        raise NotImplementedError

    def _read_particle_data_file(self, data_file, ptf, selector=None):
        # each frontend needs to implement this: read from a data_file object
        # and return a dict of fields for that data_file
        raise NotImplementedError

    def _read_particle_selection(
        self, chunks, selector, fields: List[Tuple[str, str]]
    ) -> Dict[Tuple[str, str], np.ndarray]:

        # re-organize the fields to read
        ptf, field_maps = self._get_particle_field_containers(fields)
        data: dict = {field: [] for field in fields}
        field_sizes = {field: 0 for field in fields}

        # Now we read.
        for field_r, vals in self._read_particle_fields(chunks, ptf, selector):
            # Note that we now need to check the mappings
            for field_f in field_maps[field_r]:
                data[field_f].append(vals)
                field_sizes[field_f] += vals.size

        # now concatenate into final dictionary
        rv: Dict[Tuple[str, str], np.ndarray] = {}  # the return dictionary
        fields = list(data.keys())
        for field_f in fields:
            if field_sizes[field_f]:
                vals = data.pop(field_f)
                rv[field_f] = np.concatenate(vals, axis=0).astype("float64")
            else:
                # We need to ensure the arrays have the right shape if there are no
                # particles in them.
                shape = [0]
                fname = field_f[1]
                if fname in self._vector_fields:
                    shape.append(self._vector_fields[fname])
                elif fname in self._array_fields:
                    shape.append(self._array_fields[fname])
                rv[field_f] = np.empty(shape, dtype="float64")
        return rv

    def _get_particle_field_containers(self, fields: List[Tuple[str, str]]):
        # returns ptf and field_maps:
        # ptf (particle field types) maps particle type to list of on-disk fields to read
        # field_maps stores fields, accounting for field unions

        # here we organize the fields to read by particle type, accounting for
        # particle unions

        ptf: DefaultDict[str, List[str]] = defaultdict(list)
        field_maps: DefaultDict[Tuple[str, str], List[Tuple[str, str]]] = defaultdict(
            list
        )

        unions = self.ds.particle_unions
        for field in fields:
            ftype, fname = field
            # We should add a check for p.fparticle_unions or something here
            if ftype in unions:
                for pt in unions[ftype]:
                    ptf[pt].append(fname)
                    field_maps[pt, fname].append(field)
            else:
                ptf[ftype].append(fname)
                field_maps[field].append(field)

        return ptf, field_maps

    def _read_particle_fields(self, chunks, ptf, selector):
        # this particle reader is dask-optional. If dask is enabled, it
        # constructs a dask-delayed graph and immediately computes.
        # If dask is not enabled, all dask operations are passthrough
        # functions that have no effect.

        dask_delay = _get_dask_delayed(self._dask_enabled)
        dask_compute = _get_dask_compute(self._dask_enabled)

        # Now we read each chunk
        delayed_chunks = []
        for data_file in self._sorted_chunk_iterator(chunks):
            # note: if dask is not enabled, dask_delay_wrapper is an empty
            # decorator and results are immediately returned
            read_f = self._read_particle_data_file
            delayed_chunk = dask_delay(read_f)(data_file, ptf, selector)
            delayed_chunks.append(delayed_chunk)

        # since we are returning in-mem arrays with this function, we immediately compute
        # if dask is not enabled, this will pass through the chunks without modification
        data_by_chunk = dask_compute(*delayed_chunks)

        # yield from our in-mem chunks
        for data_file_data in data_by_chunk:
            yield from data_file_data.items()

    def _sorted_chunk_iterator(self, chunks):
        chunks = list(chunks)
        data_files = set()
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        yield from sorted(data_files, key=lambda x: (x.filename, x.start))


# As a note: we don't *actually* want this to be how it is forever.  There's no
# reason we need to have the fluid and particle IO handlers separated.  But,
# for keeping track of which frontend is which, this is a useful abstraction.
class BaseParticleIOHandler(BaseIOHandler):
    pass


class IOHandlerExtracted(BaseIOHandler):

    _dataset_type = "extracted"

    def _read_data_set(self, grid, field):
        return grid.base_grid[field] / grid.base_grid.convert(field)

    def _read_data_slice(self, grid, field, axis, coord):
        sl = [slice(None), slice(None), slice(None)]
        sl[axis] = slice(coord, coord + 1)
        return grid.base_grid[field][tuple(sl)] / grid.base_grid.convert(field)
