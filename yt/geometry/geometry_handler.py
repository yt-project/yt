"""
Geometry container base class.




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import os
import cPickle
import weakref
import h5py
from exceptions import IOError, TypeError
from types import ClassType
import numpy as np
import abc
import copy

from yt.funcs import *
from yt.config import ytcfg
from yt.data_objects.field_info_container import \
    NullFunc
from yt.data_objects.particle_fields import \
    particle_deposition_functions
from yt.utilities.io_handler import io_registry
from yt.utilities.logger import ytLogger as mylog
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    ParallelAnalysisInterface, parallel_splitter
from yt.utilities.exceptions import YTFieldNotFound

def _unsupported_object(pf, obj_name):
    def _raise_unsupp(*args, **kwargs):
        raise YTObjectNotImplemented(pf, obj_name)
    return _raise_unsupp

class Index(ParallelAnalysisInterface):
    _global_mesh = True

    def __init__(self, pf, dataset_type):
        ParallelAnalysisInterface.__init__(self)
        self.parameter_file = weakref.proxy(pf)
        self.pf = self.parameter_file

        self._initialize_state_variables()

        mylog.debug("Initializing geometry index.")
        self._setup_geometry()

    def _initialize_state_variables(self):
        # TODO: Remove?
        self._max_locations = {}
        self.num_grids = None

    def _chunk(self, dobj, chunking_style, ngz = 0, **kwargs):
        # A chunk is either None or (grids, size)
        if dobj._current_chunk is None:
            self._identify_base_chunk(dobj)
        if ngz != 0 and chunking_style != "spatial":
            raise NotImplementedError
        if chunking_style == "all":
            return self._chunk_all(dobj, **kwargs)
        elif chunking_style == "spatial":
            return self._chunk_spatial(dobj, ngz, **kwargs)
        elif chunking_style == "io":
            return self._chunk_io(dobj, **kwargs)
        else:
            raise NotImplementedError

def cached_property(func):
    n = '_%s' % func.func_name
    def cached_func(self):
        if self._cache and getattr(self, n, None) is not None:
            return getattr(self, n)
        if self.data_size is None:
            tr = self._accumulate_values(n[1:])
        else:
            tr = func(self)
        if self._cache:
            setattr(self, n, tr)
        return tr
    return property(cached_func)

class YTDataChunk(object):

    def __init__(self, dobj, chunk_type, objs, data_size = None,
                 field_type = None, cache = False):
        self.dobj = dobj
        self.chunk_type = chunk_type
        self.objs = objs
        self.data_size = data_size
        self._field_type = field_type
        self._cache = cache

    def _accumulate_values(self, method):
        # We call this generically.  It's somewhat slower, since we're doing
        # costly getattr functions, but this allows us to generalize.
        mname = "select_%s" % method
        arrs = []
        for obj in self.objs:
            f = getattr(obj, mname)
            arrs.append(f(self.dobj))
        arrs = np.concatenate(arrs)
        self.data_size = arrs.shape[0]
        return arrs

    @cached_property
    def fcoords(self):
        ci = np.empty((self.data_size, 3), dtype='float64')
        if self.data_size == 0: return ci
        ind = 0
        for obj in self.objs:
            c = obj.select_fcoords(self.dobj)
            if c.shape[0] == 0: continue
            ci[ind:ind+c.shape[0], :] = c
            ind += c.shape[0]
        return ci

    @cached_property
    def icoords(self):
        ci = np.empty((self.data_size, 3), dtype='int64')
        if self.data_size == 0: return ci
        ind = 0
        for obj in self.objs:
            c = obj.select_icoords(self.dobj)
            if c.shape[0] == 0: continue
            ci[ind:ind+c.shape[0], :] = c
            ind += c.shape[0]
        return ci

    @cached_property
    def fwidth(self):
        ci = np.empty((self.data_size, 3), dtype='float64')
        if self.data_size == 0: return ci
        ind = 0
        for obj in self.objs:
            c = obj.select_fwidth(self.dobj)
            if c.shape[0] == 0: continue
            ci[ind:ind+c.shape[0], :] = c
            ind += c.shape[0]
        return ci

    @cached_property
    def ires(self):
        ci = np.empty(self.data_size, dtype='int64')
        if self.data_size == 0: return ci
        ind = 0
        for obj in self.objs:
            c = obj.select_ires(self.dobj)
            if c.shape == 0: continue
            ci[ind:ind+c.size] = c
            ind += c.size
        return ci

    @cached_property
    def tcoords(self):
        self.dtcoords
        return self._tcoords

    @cached_property
    def dtcoords(self):
        ct = np.empty(self.data_size, dtype='float64')
        cdt = np.empty(self.data_size, dtype='float64')
        self._tcoords = ct # Se this for tcoords
        if self.data_size == 0: return cdt
        ind = 0
        for obj in self.objs:
            gdt, gt = obj.tcoords(self.dobj)
            if gt.shape == 0: continue
            ct[ind:ind+gt.size] = gt
            cdt[ind:ind+gdt.size] = gdt
            ind += gt.size
        return cdt

class ChunkDataCache(object):
    def __init__(self, base_iter, preload_fields, geometry_handler,
                 max_length = 256):
        # At some point, max_length should instead become a heuristic function,
        # potentially looking at estimated memory usage.  Note that this never
        # initializes the iterator; it assumes the iterator is already created,
        # and it calls next() on it.
        self.base_iter = base_iter.__iter__()
        self.queue = []
        self.max_length = max_length
        self.preload_fields = preload_fields
        self.geometry_handler = geometry_handler
        self.cache = {}

    def __iter__(self):
        return self
    
    def next(self):
        if len(self.queue) == 0:
            for i in range(self.max_length):
                try:
                    self.queue.append(self.base_iter.next())
                except StopIteration:
                    break
            # If it's still zero ...
            if len(self.queue) == 0: raise StopIteration
            chunk = YTDataChunk(None, "cache", self.queue, cache=False)
            self.cache = self.geometry_handler.io._read_chunk_data(
                chunk, self.preload_fields)
        g = self.queue.pop(0)
        g._initialize_cache(self.cache.pop(g.id, {}))
        return g
