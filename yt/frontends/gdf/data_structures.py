"""
Data structures for Chombo.

Author: Matthew Turk <matthewturk@gmail.com>
Author: J. S. Oishi <jsoishi@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/
License:
  Copyright (C) 2008-2011 Matthew Turk, J. S. Oishi.  All Rights Reserved.

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

from yt.funcs import *
from yt.data_objects.grid_patch import \
           AMRGridPatch
from yt.data_objects.hierarchy import \
           AMRHierarchy
from yt.data_objects.static_output import \
           StaticOutput

from .fields import GDFFieldContainer

class GDFGrid(AMRGridPatch):
    _id_offset = 0
    def __init__(self, id, hierarchy, level, start, dimensions):
        AMRGridPatch.__init__(self, id, filename = hierarchy.hierarchy_filename,
                              hierarchy = hierarchy)
        self.Parent = []
        self.Children = []
        self.Level = level
        self.start_index = start.copy()
        self.stop_index = self.start_index + dimensions
        self.ActiveDimensions = dimensions.copy()

    def _setup_dx(self):
        # So first we figure out what the index is.  We don't assume
        # that dx=dy=dz , at least here.  We probably do elsewhere.
        id = self.id - self._id_offset
        if len(self.Parent) > 0:
            self.dds = self.Parent[0].dds / self.pf.refine_by
        else:
            LE, RE = self.hierarchy.grid_left_edge[id,:], \
                     self.hierarchy.grid_right_edge[id,:]
            self.dds = na.array((RE-LE)/self.ActiveDimensions)
        if self.pf.dimensionality < 2: self.dds[1] = 1.0
        if self.pf.dimensionality < 3: self.dds[2] = 1.0
        self.data['dx'], self.data['dy'], self.data['dz'] = self.dds

class GDFHierarchy(AMRHierarchy):

    grid = GDFGrid
    
    def __init__(self, pf, data_style='grid_data_format'):
        self.parameter_file = weakref.proxy(pf)
        # for now, the hierarchy file is the parameter file!
        self.hierarchy_filename = self.parameter_file.parameter_filename
        self.directory = os.path.dirname(self.hierarchy_filename)
        self._fhandle = h5py.File(self.hierarchy_filename)
        AMRHierarchy.__init__(self,pf,data_style)

        self._fhandle.close()

    def _initialize_data_storage(self):
        pass

    def _detect_fields(self):
        ncomp = int(self._fhandle['/'].attrs['num_components'])
        self.field_list = [c[1] for c in self._fhandle['/'].attrs.listitems()[-ncomp:]]
    
    def _setup_classes(self):
        dd = self._get_data_reader_dict()
        AMRHierarchy._setup_classes(self, dd)
        self.object_types.sort()

    def _count_grids(self):
        self.num_grids = 0
        for lev in self._levels:
            self.num_grids += self._fhandle[lev]['Processors'].len()
        
    def _parse_hierarchy(self):
        f = self._fhandle # shortcut
        
        # this relies on the first Group in the H5 file being
        # 'Chombo_global'
        levels = f.listnames()[1:]
        self.grids = []
        i = 0
        for lev in levels:
            level_number = int(re.match('level_(\d+)',lev).groups()[0])
            boxes = f[lev]['boxes'].value
            dx = f[lev].attrs['dx']
            for level_id, box in enumerate(boxes):
                si = na.array([box['lo_%s' % ax] for ax in 'ijk'])
                ei = na.array([box['hi_%s' % ax] for ax in 'ijk'])
                pg = self.grid(len(self.grids),self,level=level_number,
                               start = si, stop = ei)
                self.grids.append(pg)
                self.grids[-1]._level_id = level_id
                self.grid_left_edge[i] = dx*si.astype(self.float_type)
                self.grid_right_edge[i] = dx*(ei.astype(self.float_type) + 1)
                self.grid_particle_count[i] = 0
                self.grid_dimensions[i] = ei - si + 1
                i += 1
        temp_grids = na.empty(len(grids), dtype='object')
        for gi, g in enumerate(self.grids): temp_grids[gi] = g
        self.grids = temp_grids

    def _populate_grid_objects(self):
        for g in self.grids:
            g._prepare_grid()
            g._setup_dx()

        for g in self.grids:
            g.Children = self._get_grid_children(g)
            for g1 in g.Children:
                g1.Parent.append(g)
        self.max_level = self.grid_levels.max()

    def _setup_unknown_fields(self):
        pass

    def _setup_derived_fields(self):
        self.derived_field_list = []

    def _get_grid_children(self, grid):
        mask = na.zeros(self.num_grids, dtype='bool')
        grids, grid_ind = self.get_box_grids(grid.LeftEdge, grid.RightEdge)
        mask[grid_ind] = True
        return [g for g in self.grids[mask] if g.Level == grid.Level + 1]

class GDFStaticOutput(StaticOutput):
    _hierarchy_class = GDFHierarchy
    _fieldinfo_class = GDFFieldContainer
    
    def __init__(self, filename, data_style='grid_data_format',
                 storage_filename = None):
        StaticOutput.__init__(self, filename, data_style)
        self._handle = h5py.File(self.filename, "r")
        self.storage_filename = storage_filename
        self.field_info = self._fieldinfo_class()
        self._handle.close()
        del self._handle
        
    def _set_units(self):
        """
        Generates the conversion to various physical _units based on the parameter file
        """
        self.units = {}
        self.time_units = {}
        if len(self.parameters) == 0:
            self._parse_parameter_file()
        self.time_units['1'] = 1
        self.units['1'] = 1.0
        self.units['unitary'] = 1.0 / (self.domain_right_edge - self.domain_right_edge).max()
        seconds = 1
        self.time_units['years'] = seconds / (365*3600*24.0)
        self.time_units['days']  = seconds / (3600*24.0)
        # This should be improved.
        for field_name in self._handle["/field_types"]:
            self.units[field_name] = self._handle["/%s/field_to_cgs" % field_name]

    def _parse_parameter_file(self):
        sp = self._handle["/simulation_parameters"].attrs
        self.domain_left_edge = sp["domain_left_edge"][:]
        self.domain_right_edge = sp["domain_right_edge"][:]
        self.refine_by = sp["refine_by"][:]
        self.dimensionality = sp["dimensionality"][:]
        self.current_time = sp["current_time"][:]
        self.unique_identifier = sp["unique_identifier"]
        self.cosmological_simulation = sp["cosmological_simulation"]
        if sp["num_ghost_zones"] != 0: raise RuntimeError
        self.field_ordering = sp["field_ordering"]
        self.boundary_conditions = sp["boundary_conditions"][:]
        if self.cosmological_simulation:
            self.current_redshift = sp["current_redshift"]
            self.omega_lambda = sp["omega_lambda"]
            self.omega_matter = sp["omega_matter"]
            self.hubble_constant = sp["hubble_constant"]
        else:
            self.current_redshift = self.omega_lambda = self.omega_matter = \
                self.hubble_constant = self.cosmological_simulation = 0.0
        
    @classmethod
    def _is_valid(self, *args, **kwargs):
        try:
            fileh = h5py.File(args[0],'r')
            if "gridded_data_format" in fileh:
                return True
        except:
            pass
        return False


