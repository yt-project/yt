"""
Data structures for Streaming, in-memory datasets

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt-project.org/
License:
  Copyright (C) 2011 Matthew Turk.  All Rights Reserved.

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

import weakref
import numpy as na

from yt.utilities.io_handler import io_registry
from yt.funcs import *
from yt.config import ytcfg
from yt.data_objects.grid_patch import \
    AMRGridPatch
from yt.geometry.grid_geometry_handler import \
    GridGeometryHandler
from yt.data_objects.static_output import \
    StaticOutput
from yt.utilities.logger import ytLogger as mylog
from yt.data_objects.field_info_container import \
    FieldInfoContainer, NullFunc
from yt.utilities.lib import \
    get_box_grids_level

from .fields import \
    StreamFieldInfo, \
    add_stream_field, \
    KnownStreamFields

class StreamGrid(AMRGridPatch):
    """
    Class representing a single In-memory Grid instance.
    """

    __slots__ = ['proc_num']
    _id_offset = 0
    def __init__(self, id, hierarchy):
        """
        Returns an instance of StreamGrid with *id*, associated with *filename*
        and *hierarchy*.
        """
        #All of the field parameters will be passed to us as needed.
        AMRGridPatch.__init__(self, id, filename = None, hierarchy = hierarchy)
        self._children_ids = []
        self._parent_id = -1
        self.Level = -1

    def _guess_properties_from_parent(self):
        rf = self.pf.refine_by
        my_ind = self.id - self._id_offset
        le = self.LeftEdge
        self.dds = self.Parent.dds/rf
        ParentLeftIndex = na.rint((self.LeftEdge-self.Parent.LeftEdge)/self.Parent.dds)
        self.start_index = rf*(ParentLeftIndex + self.Parent.get_global_startindex()).astype('int64')
        self.LeftEdge = self.Parent.LeftEdge + self.Parent.dds * ParentLeftIndex
        self.RightEdge = self.LeftEdge + self.ActiveDimensions*self.dds
        self.hierarchy.grid_left_edge[my_ind,:] = self.LeftEdge
        self.hierarchy.grid_right_edge[my_ind,:] = self.RightEdge
        self._child_mask = None
        self._child_index_mask = None
        self._child_indices = None
        self._setup_dx()

    def set_filename(self, filename):
        pass

    def __repr__(self):
        return "StreamGrid_%04i" % (self.id)

    @property
    def Parent(self):
        if self._parent_id == -1: return None
        return self.hierarchy.grids[self._parent_id - self._id_offset]

    @property
    def Children(self):
        return [self.hierarchy.grids[cid - self._id_offset]
                for cid in self._children_ids]

class StreamHandler(object):
    def __init__(self, left_edges, right_edges, dimensions,
                 levels, parent_ids, particle_count, processor_ids,
                 fields, io = None):
        self.left_edges = left_edges
        self.right_edges = right_edges
        self.dimensions = dimensions
        self.levels = levels
        self.parent_ids = parent_ids
        self.particle_count = particle_count
        self.processor_ids = processor_ids
        self.num_grids = self.levels.size
        self.fields = fields
        self.io = io

    def get_fields(self):
        return self.fields.all_fields

class StreamHierarchy(GridGeometryHandler):

    grid = StreamGrid

    def __init__(self, pf, data_style = None):
        self.data_style = data_style
        self.float_type = 'float64'
        self.parameter_file = weakref.proxy(pf) # for _obtain_enzo
        self.stream_handler = pf.stream_handler
        self.float_type = "float64"
        self.directory = os.getcwd()
        GridGeometryHandler.__init__(self, pf, data_style)

    def _initialize_data_storage(self):
        pass

    def _count_grids(self):
        self.num_grids = self.stream_handler.num_grids

    def _setup_unknown_fields(self):
        for field in self.field_list:
            if field in self.parameter_file.field_info: continue
            mylog.info("Adding %s to list of fields", field)
            cf = None
            if self.parameter_file.has_key(field):
                def external_wrapper(f):
                    def _convert_function(data):
                        return data.convert(f)
                    return _convert_function
                cf = external_wrapper(field)
            # Note that we call add_field on the field_info directly.  This
            # will allow the same field detection mechanism to work for 1D, 2D
            # and 3D fields.
            self.pf.field_info.add_field(
                    field, lambda a, b: None,
                    convert_function=cf, take_log=False)
            

    def _parse_hierarchy(self):
        self.grid_dimensions = self.stream_handler.dimensions
        self.grid_left_edge[:] = self.stream_handler.left_edges
        self.grid_right_edge[:] = self.stream_handler.right_edges
        self.grid_levels[:] = self.stream_handler.levels
        self.grid_procs = self.stream_handler.processor_ids
        self.grid_particle_count[:] = self.stream_handler.particle_count
        mylog.debug("Copying reverse tree")
        self.grids = []
        # We enumerate, so it's 0-indexed id and 1-indexed pid
        for id in xrange(self.num_grids):
            self.grids.append(self.grid(id, self))
            self.grids[id].Level = self.grid_levels[id, 0]
        parent_ids = self.stream_handler.parent_ids
        if parent_ids is not None:
            reverse_tree = self.stream_handler.parent_ids.tolist()
            # Initial setup:
            for gid,pid in enumerate(reverse_tree):
                if pid >= 0:
                    self.grids[id]._parent_id = pid
                    self.grids[pid]._children_ids.append(self.grids[gid].id)
        else:
            mylog.debug("Reconstructing parent-child relationships")
            self._reconstruct_parent_child()
        self.max_level = self.grid_levels.max()
        mylog.debug("Preparing grids")
        temp_grids = na.empty(self.num_grids, dtype='object')
        for i, grid in enumerate(self.grids):
            if (i%1e4) == 0: mylog.debug("Prepared % 7i / % 7i grids", i, self.num_grids)
            grid.filename = None
            grid._prepare_grid()
            grid.proc_num = self.grid_procs[i]
            temp_grids[i] = grid
        self.grids = temp_grids
        mylog.debug("Prepared")

    def _reconstruct_parent_child(self):
        mask = na.empty(len(self.grids), dtype='int32')
        mylog.debug("First pass; identifying child grids")
        for i, grid in enumerate(self.grids):
            get_box_grids_level(self.grid_left_edge[i,:],
                                self.grid_right_edge[i,:],
                                self.grid_levels[i] + 1,
                                self.grid_left_edge, self.grid_right_edge,
                                self.grid_levels, mask)
            ids = na.where(mask.astype("bool"))
            grid._children_ids = ids[0] # where is a tuple
        mylog.debug("Second pass; identifying parents")
        for i, grid in enumerate(self.grids): # Second pass
            for child in grid.Children:
                child._parent_id = i

    def _initialize_grid_arrays(self):
        GridGeometryHandler._initialize_grid_arrays(self)
        self.grid_procs = na.zeros((self.num_grids,1),'int32')

    def save_data(self, *args, **kwargs):
        pass

    def _detect_fields(self):
        self.field_list = list(set(self.stream_handler.get_fields()))

    def _setup_derived_fields(self):
        self.derived_field_list = []
        for field in self.parameter_file.field_info:
            try:
                fd = self.parameter_file.field_info[field].get_dependencies(
                            pf = self.parameter_file)
            except:
                continue
            available = na.all([f in self.field_list for f in fd.requested])
            if available: self.derived_field_list.append(field)
        for field in self.field_list:
            if field not in self.derived_field_list:
                self.derived_field_list.append(field)

    def _setup_classes(self):
        dd = self._get_data_reader_dict()
        GridGeometryHandler._setup_classes(self, dd)
        self.object_types.sort()

    def _populate_grid_objects(self):
        for g in self.grids:
            g._setup_dx()
        self.max_level = self.grid_levels.max()

    def _setup_data_io(self):
        if self.stream_handler.io is not None:
            self.io = self.stream_handler.io
        else:
            self.io = io_registry[self.data_style](self.stream_handler)

class StreamStaticOutput(StaticOutput):
    _hierarchy_class = StreamHierarchy
    _fieldinfo_fallback = StreamFieldInfo
    _fieldinfo_known = KnownStreamFields
    _data_style = 'stream'

    def __init__(self, stream_handler):
        #if parameter_override is None: parameter_override = {}
        #self._parameter_override = parameter_override
        #if conversion_override is None: conversion_override = {}
        #self._conversion_override = conversion_override

        self.stream_handler = stream_handler
        StaticOutput.__init__(self, "InMemoryParameterFile", self._data_style)

        self.units = {}
        self.time_units = {}

    def _parse_parameter_file(self):
        self.basename = self.stream_handler.name
        self.parameters['CurrentTimeIdentifier'] = time.time()
        self.unique_identifier = self.parameters["CurrentTimeIdentifier"]
        self.domain_left_edge = self.stream_handler.domain_left_edge[:]
        self.domain_right_edge = self.stream_handler.domain_right_edge[:]
        self.refine_by = self.stream_handler.refine_by
        self.dimensionality = self.stream_handler.dimensionality
        self.domain_dimensions = self.stream_handler.domain_dimensions
        self.current_time = self.stream_handler.simulation_time
        if self.stream_handler.cosmology_simulation:
            self.cosmological_simulation = 1
            self.current_redshift = self.stream_handler.current_redshift
            self.omega_lambda = self.stream_handler.omega_lambda
            self.omega_matter = self.stream_handler.omega_matter
            self.hubble_constant = self.stream_handler.hubble_constant
        else:
            self.current_redshift = self.omega_lambda = self.omega_matter = \
                self.hubble_constant = self.cosmological_simulation = 0.0

    def _set_units(self):
        pass

    @classmethod
    def _is_valid(cls, *args, **kwargs):
        return False
