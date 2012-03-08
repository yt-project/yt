"""
Handler for Octree-based geometries

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/
License:
  Copyright (C) 2007-2011 Matthew Turk.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import h5py
import numpy as na
import string, re, gc, time, cPickle
import weakref

from itertools import chain, izip

from yt.funcs import *
from yt.utilities.logger import ytLogger as mylog
from yt.arraytypes import blankRecordArray
from yt.config import ytcfg
from yt.data_objects.field_info_container import NullFunc
from yt.geometry.geometry_handler import GeometryHandler, YTDataChunk
from yt.utilities.definitions import MAXLEVEL
from yt.utilities.io_handler import io_registry
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    ParallelAnalysisInterface, parallel_splitter
from object_finding_mixin import ObjectFindingMixin

from yt.data_objects.data_containers import data_object_registry

class GridGeometryHandler(ObjectFindingMixin, GeometryHandler):
    float_type = 'float64'

    def _setup_geometry(self):
        mylog.debug("Counting octs.")
        self._count_octs()

        mylog.debug("Initializing Octree structure.")
        self._initialize_octree()

        mylog.debug("Re-examining hierarchy")
        self._initialize_stats()

    def _create_octree_structure(self):
        self.oct_levels = na.zeros((self.num_octs,1), dtype='int32')
        self.oct_file_locations = na.zeros((self.num_octs, self._nfl), dtype='int64')
        self.oct_hilbert_indices = na.zeros(self.num_octs, dtype='uint64')
        # 8 bits, 8 octs:
        self.oct_refined = na.zeros((self.num_octs, 1), dtype='int8')
    
    def convert(self, unit):
        return self.parameter_file.conversion_factors[unit]

    def _identify_base_chunk(self, dobj):
        if getattr(dobj, "_octs", None) is None:
            oct_mask = dobj.selector.select_octs(self.max_level, self.oct_hilbert_indices)
            dobj._octs = oct_mask
        if getattr(dobj, "size", None) is None:
            dobj.size = self._count_selection(dobj)
            dobj.shape = (dobj.size,)
        dobj._current_chunk = list(self._chunk_all(dobj))[0]

    def _count_selection(self, dobj, oct_mask = None):
        if oct_mask is None: oct_mask = dobj._octs

        count = ((self.count(dobj.selector) for g in grids))
        return count

    def _read_particle_fields(self, fields, dobj, chunk = None):
        if len(fields) == 0: return {}, []
        selector = dobj.selector
        if chunk is None:
            self._identify_base_chunk(dobj)
        fields_to_return = {}
        fields_to_read, fields_to_generate = [], []
        for ftype, fname in fields:
            if fname in self.field_list:
                fields_to_read.append((ftype, fname))
            else:
                fields_to_generate.append((ftype, fname))
        if len(fields_to_read) == 0:
            return {}, fields_to_generate
        fields_to_return = self.io._read_particle_selection(
                    self._chunk_io(dobj), selector,
                    fields_to_read)
        for field in fields_to_read:
            ftype, fname = field
            conv_factor = self.pf.field_info[fname]._convert_function(self)
            na.multiply(fields_to_return[field], conv_factor,
                        fields_to_return[field])
        return fields_to_return, fields_to_generate

    def _read_fluid_fields(self, fields, dobj, chunk = None):
        if len(fields) == 0: return {}, []
        selector = dobj.selector
        if chunk is None:
            self._identify_base_chunk(dobj)
            chunk_size = dobj.size
        else:
            chunk_size = chunk.data_size
        fields_to_return = {}
        fields_to_read, fields_to_generate = [], []
        for ftype, fname in fields:
            if fname in self.field_list:
                fields_to_read.append((ftype, fname))
            else:
                fields_to_generate.append((ftype, fname))
        if len(fields_to_read) == 0:
            return {}, fields_to_generate
        fields_to_return = self.io._read_fluid_selection(self._chunk_io(dobj),
                                                   selector,
                                                   fields_to_read,
                                                   chunk_size)
        for field in fields_to_read:
            ftype, fname = field
            conv_factor = self.pf.field_info[fname]._convert_function(self)
            na.multiply(fields_to_return[field], conv_factor,
                        fields_to_return[field])
        #mylog.debug("Don't know how to read %s", fields_to_generate)
        return fields_to_return, fields_to_generate

    def _chunk_all(self, dobj):
        gobjs = getattr(dobj._current_chunk, "objs", dobj._grids)
        yield YTDataChunk(dobj, "all", gobjs, dobj.size)
        
    def _chunk_spatial(self, dobj, ngz):
        gobjs = getattr(dobj._current_chunk, "objs", dobj._grids)
        for i,og in enumerate(gobjs):
            if ngz > 0:
                g = og.retrieve_ghost_zones(ngz, [], smoothed=True)
            else:
                g = og
            size = self._count_selection(dobj, [og])
            if size == 0: continue
            yield YTDataChunk(dobj, "spatial", [g], size)

    def _chunk_io(self, dobj):
        gfiles = defaultdict(list)
        gobjs = getattr(dobj._current_chunk, "objs", dobj._grids)
        for g in gobjs:
            gfiles[g.filename].append(g)
        for fn in sorted(gfiles):
            gs = gfiles[fn]
            yield YTDataChunk(dobj, "io", gs, self._count_selection(dobj, gs))
