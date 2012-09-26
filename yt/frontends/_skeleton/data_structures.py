"""
Skeleton data structures

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt-project.org/
License:
  Copyright (C) 2012 Matthew Turk.  All Rights Reserved.

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
import stat
import numpy as np
import weakref

from yt.funcs import *
from yt.data_objects.grid_patch import \
    AMRGridPatch
from yt.data_objects.hierarchy import \
    AMRHierarchy
from yt.data_objects.static_output import \
    StaticOutput
from yt.utilities.definitions import \
    mpc_conversion, sec_conversion
from yt.utilities.io_handler import \
    io_registry
from yt.utilities.physical_constants import cm_per_mpc
from .fields import SkeletonFieldInfo, add_flash_field, KnownSkeletonFields
from yt.data_objects.field_info_container import \
    FieldInfoContainer, NullFunc, ValidateDataField, TranslationFunc

class SkeletonGrid(AMRGridPatch):
    _id_offset = 0
    #__slots__ = ["_level_id", "stop_index"]
    def __init__(self, id, hierarchy, level):
        AMRGridPatch.__init__(self, id, filename = hierarchy.hierarchy_filename,
                              hierarchy = hierarchy)
        self.Parent = None
        self.Children = []
        self.Level = level

    def __repr__(self):
        return "SkeletonGrid_%04i (%s)" % (self.id, self.ActiveDimensions)

class SkeletonHierarchy(AMRHierarchy):

    grid = SkeletonGrid
    float_type = np.float64
    
    def __init__(self, pf, data_style='skeleton'):
        self.data_style = data_style
        self.parameter_file = weakref.proxy(pf)
        # for now, the hierarchy file is the parameter file!
        self.hierarchy_filename = self.parameter_file.parameter_filename
        self.directory = os.path.dirname(self.hierarchy_filename)
        AMRHierarchy.__init__(self, pf, data_style)

    def _initialize_data_storage(self):
        pass

    def _detect_fields(self):
        # This needs to set a self.field_list that contains all the available,
        # on-disk fields.
        pass
    
    def _count_grids(self):
        # This needs to set self.num_grids
        pass
        
    def _parse_hierarchy(self):
        # This needs to fill the following arrays, where N is self.num_grids:
        #   self.grid_left_edge         (N, 3) <= float64
        #   self.grid_right_edge        (N, 3) <= float64
        #   self.grid_dimensions        (N, 3) <= int
        #   self.grid_particle_count    (N, 1) <= int
        #   self.grid_levels            (N, 1) <= int
        #   self.grids                  (N, 1) <= grid objects
        #   
        pass
                        
    def _populate_grid_objects(self):
        # For each grid, this must call:
        #   grid._prepare_grid()
        #   grid._setup_dx()
        # This must also set:
        #   grid.Children <= list of child grids
        #   grid.Parent   <= parent grid
        # This is handled by the frontend because often the children must be
        # identified.
        pass

class SkeletonStaticOutput(StaticOutput):
    _hierarchy_class = SkeletonHierarchy
    _fieldinfo_fallback = SkeletonFieldInfo
    _fieldinfo_known = KnownSkeletonFields
    _handle = None
    
    def __init__(self, filename, data_style='skeleton',
                 storage_filename = None,
                 conversion_override = None):

        if conversion_override is None: conversion_override = {}
        self._conversion_override = conversion_override

        StaticOutput.__init__(self, filename, data_style)
        self.storage_filename = storage_filename

    def _set_units(self):
        # This needs to set up the dictionaries that convert from code units to
        # CGS.  The needed items are listed in the second entry:
        #   self.time_units         <= sec_conversion
        #   self.conversion_factors <= mpc_conversion
        #   self.units              <= On-disk fields
        pass

    def _parse_parameter_file(self):
        # This needs to set up the following items:
        #
        #   self.unique_identifier
        #   self.parameters             <= full of code-specific items of use
        #   self.domain_left_edge       <= array of float64
        #   self.domain_right_edge      <= array of float64
        #   self.dimensionality         <= int
        #   self.domain_dimensions      <= array of int64
        #   self.current_time           <= simulation time in code units
        #
        # We also set up cosmological information.  Set these to zero if
        # non-cosmological.
        #
        #   self.cosmological_simulation    <= int, 0 or 1
        #   self.current_redshift           <= float
        #   self.omega_lambda               <= float
        #   self.omega_matter               <= float
        #   self.hubble_constant            <= float

    @classmethod
    def _is_valid(self, *args, **kwargs):
        # This accepts a filename or a set of arguments and returns True or
        # False depending on if the file is of the type requested.
        return False


