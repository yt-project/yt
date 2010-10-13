"""
ART-specific data structures

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: UCSD
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2010 Matthew Turk.  All Rights Reserved.

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

import numpy as na
import stat
import weakref
import cPickle
import art_reader

from yt.funcs import *
from yt.data_objects.grid_patch import \
      AMRGridPatch
from yt.data_objects.hierarchy import \
      AMRHierarchy
from yt.data_objects.static_output import \
      StaticOutput
from .fields import ARTFieldContainer
from yt.utilities.definitions import \
    mpc_conversion
from yt.utilities.io_handler import \
    io_registry

class ARTGrid(AMRGridPatch):
    _id_offset = 0
    #__slots__ = ["_level_id", "stop_index"]
    def __init__(self, id, hierarchy, LE,RE,Dim,Level,Parent ):
        AMRGridPatch.__init__(self, id, filename = hierarchy.hierarchy_filename,
                              hierarchy = hierarchy)
        self.id = id
        self.LE,self.RE = LE,RE
        self.Level = Level
        self.Parent = Parent
        self.Children = []

    def __repr__(self):
        return "ARTGrid_%04i (%s)" % (self.id, self.ActiveDimensions)

class ARTHierarchy(AMRHierarchy):

    grid = ARTGrid
    _handle = None
    
    def __init__(self,pf,data_style='art'):
        AMRHierarchy.__init__(self,pf,data_style)

    def _detect_fields(self):
        self.field_list = ['gas_density','total_density','px','py','pz'
            'energy','pressure','gamma','internale','pot_new','pot_old']
        
    def _setup_classes(self):
        dd = self._get_data_reader_dict()
        AMRHierarchy._setup_classes(self, dd)
        self.object_types.sort()
        
    def _count_grids(self):
        self.num_grids = len(self.pf.art.grid_level)
        #the length of any grid array should suffice
        
    def _parse_hierarchy(self):
        #all of the grid edges=indices are defined on the finest Level
        #what ART calls the root grid. To transform indices and dimensionality
        #to the relative level of a cell/grid, divide by the dimensions of the 
        # grid. Because ART is totally octree, and every cell on a level is the 
        #same size, there are no 'colaseced' regions 
        art = self.pf.art #alias
        
        self.grid_left_edge  = art.grid_left_index / art.grid_dimensions
        self.grid_right_edge = (art.grid_left_index + \
                               art.grid_dimensions) / art.grid_dimensions
        self.grid_dimensions = art.grid_dimensions
        self.grid_levels = art.grid_level
        self.grid_particle_count = art.grid_level*0
        self.grid_parents = art.grid_parents
        
        
    def _populate_grid_objects(self):
        iters = zip(range(self.num_grids),self.grid_left_edge,
                self.grid_right_edge,self.grid_dimensions,self.grid_levels,
                self.grid_parents)
        grids = []
        for idx,LE,RE,Dim,Level,Parent in iters:
            g = self.grid(idx,self,LE,RE,Dim,Level,Parent)
            g._prepare_grid()
            g._setup_dx()
            grids.append(g)
        self.grids = na.array(grids,dtype='object')
      
    def _setup_unknown_fields(self):
        pass

    def _setup_derived_fields(self):
        pass
        
class ARTStaticOutput(StaticOutput):
    _hierarchy_class = ARTHierarchy
    _fieldinfo_class = ARTFieldContainer
    _handle = None
    
    def __init__(self, filename, data_style='art',
                 storage_filename = None):
        from art_reader import *                 
        StaticOutput.__init__(self, filename, data_style)
        self.storage_filename = storage_filename
        self.art = cPickle.load(open(filename,'rb'))
        #This is Grid class from the art_reader file
        self.field_info = self._fieldinfo_class()
        # hardcoded for now
        
        self.current_time = 0.0
        # These should be explicitly obtained from the file, but for now that
        # will wait until a reorganization of the source tree and better
        # generalization.
        self.dimensionality = 3
        self.refine_by = 2
        self.parameters["HydroMethod"] = 'art'
        self.parameters["Time"] = 1. # default unit is 1...

    def __repr__(self):
        return self.basename.rsplit(".", 1)[0]
        
    def _set_units(self):
        """
        Generates the conversion to various physical _units based on the parameter file
        """
        self.units = {}
        self.time_units = {}
        if len(self.parameters) == 0:
            self._parse_parameter_file()
        self._setup_nounits_units()
        self.conversion_factors = defaultdict(lambda: 1.0)
        self.time_units['1'] = 1
        self.units['1'] = 1.0
        self.units['unitary'] = 1.0 / (self.domain_right_edge - self.domain_left_edge).max()
        seconds = 1 #self["Time"]
        self.time_units['years'] = seconds / (365*3600*24.0)
        self.time_units['days']  = seconds / (3600*24.0)

    def _setup_nounits_units(self):
        z = 0
        mylog.warning("Setting 1.0 in code units to be 1.0 cm")
        if not self.has_key("TimeUnits"):
            mylog.warning("No time units.  Setting 1.0 = 1 second.")
            self.conversion_factors["Time"] = 1.0
        for unit in mpc_conversion.keys():
            self.units[unit] = mpc_conversion[unit] / mpc_conversion["cm"]

    @classmethod
    def _is_valid(self, *args, **kwargs):
        return os.path.exists(fn)

