"""
Data structures for Gadget.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: UCSD
Author: Chris Moody <cemoody@ucsc.edu>
Affiliation: UCSC
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
from itertools import izip

from yt.funcs import *
from yt.data_objects.grid_patch import \
    AMRGridPatch
from yt.data_objects.hierarchy import \
    AMRHierarchy
from yt.data_objects.static_output import \
    StaticOutput

from .fields import GadgetFieldInfo, KnownGadgetFields
from yt.data_objects.field_info_container import \
    FieldInfoContainer, NullFunc

class GadgetGrid(AMRGridPatch):
    _id_offset = 0
    def __init__(self, hierarchy, id, dimensions, start,
                 level, parent_id, particle_count):
        AMRGridPatch.__init__(self, id, filename = hierarchy.filename,
                              hierarchy = hierarchy)
        self.Parent = [] # Only one parent per grid        
        self.Children = []
        self.Level = level
        self.ActiveDimensions = dimensions.copy()
        self.NumberOfParticles = particle_count
        self.start_index = start.copy().astype("int64")
        self.stop_index = self.start_index + dimensions.copy()
        self.id = id
        self._parent_id = parent_id
        
        try:
            padd = '/data/grid_%010i/particles' % id
            self.particle_types = self.hierarchy._handle[padd].keys()
        except:
            self.particle_types =  ()
        self.filename = hierarchy.filename
        
    def __repr__(self):
        return 'GadgetGrid_%05i'%self.id
        
class GadgetHierarchy(AMRHierarchy):
    grid = GadgetGrid

    def __init__(self, pf, data_style='gadget_hdf5'):
        self.filename = pf.filename
        self.directory = os.path.dirname(pf.filename)
        self.data_style = data_style
        self._handle = h5py.File(pf.filename)
        AMRHierarchy.__init__(self, pf, data_style)
        self._handle.close()
        self._handle = None
        

    def _initialize_data_storage(self):
        pass

    def _detect_fields(self):
        #this adds all the fields in 
        #/particle_types/{Gas/Stars/etc.}/{position_x, etc.}
        self.field_list = []
        for ptype in self._handle['particle_types'].keys():
            for field in self._handle['particle_types'+'/'+ptype]:
                if field not in self.field_list:
                    self.field_list += field,
        
    def _setup_classes(self):
        dd = self._get_data_reader_dict()
        AMRHierarchy._setup_classes(self, dd)
        self.object_types.sort()

    def _count_grids(self):
        self.num_grids = len(self._handle['/grid_dimensions'])
        
    def _parse_hierarchy(self):
        f = self._handle # shortcut
        npa = na.array
        DLE = self.parameter_file.domain_left_edge
        DRE = self.parameter_file.domain_right_edge
        DW = (DRE - DLE)
        
        self.grid_levels.flat[:] = f['/grid_level'][:].astype("int32")
        LI = f['/grid_left_index'][:]
        print LI
        self.grid_dimensions[:] = f['/grid_dimensions'][:]
        self.grid_left_edge[:]  = (LI * DW + DLE)
        dxs = 1.0/(2**(self.grid_levels+1)) * DW
        self.grid_right_edge[:] = self.grid_left_edge \
                                + dxs *(1 + self.grid_dimensions)
        self.grid_particle_count.flat[:] = f['/grid_particle_count'][:].astype("int32")
        grid_parent_id = f['/grid_parent_id'][:]
        self.max_level = na.max(self.grid_levels)
        
        args = izip(xrange(self.num_grids), self.grid_levels.flat,
                    grid_parent_id, LI,
                    self.grid_dimensions, self.grid_particle_count.flat)
        self.grids = na.empty(len(args), dtype='object')
        for gi, (j,lvl,p, le, d, n) in enumerate(args):
            self.grids[gi] = self.grid(self,j,d,le,lvl,p,n)
        
    def _populate_grid_objects(self):    
        for g in self.grids:
            if g._parent_id >= 0:
                parent = self.grids[g._parent_id]
                g.Parent = parent
                parent.Children.append(g)
        for g in self.grids:
            g._prepare_grid()
            g._setup_dx()
            
    def _setup_derived_fields(self):
        self.derived_field_list = []

class GadgetStaticOutput(StaticOutput):
    _hierarchy_class = GadgetHierarchy
    _fieldinfo_fallback = GadgetFieldInfo
    _fieldinfo_known = KnownGadgetFields

    def __init__(self, filename,storage_filename=None) :
        self.storage_filename = storage_filename
        self.filename = filename
        
        StaticOutput.__init__(self, filename, 'gadget_infrastructure')
        self._set_units()
        
    def _set_units(self):
        self.units = {}
        self.time_units = {}
        self.time_units['1'] = 1
        self.units['1'] = 1.0
        self.units['cm'] = 1.0
        self.units['unitary'] = 1.0 / \
            (self.domain_right_edge - self.domain_left_edge).max()
        seconds = 1 #self["Time"]
        self.time_units['years'] = seconds / (365*3600*24.0)
        self.time_units['days']  = seconds / (3600*24.0)
        self.time_units['Myr'] = self.time_units['years'] / 1.0e6
        self.time_units['Gyr']  = self.time_units['years'] / 1.0e9

    def _parse_parameter_file(self):
        fileh = h5py.File(self.filename)
        sim_param = fileh['/simulation_parameters'].attrs
        self.refine_by = sim_param['refine_by']
        self.dimensionality = sim_param['dimensionality']
        self.num_ghost_zones = sim_param['num_ghost_zones']
        self.field_ordering = sim_param['field_ordering']
        self.domain_dimensions = sim_param['domain_dimensions']
        self.current_time = sim_param['current_time']
        self.domain_left_edge = sim_param['domain_left_edge']
        self.domain_right_edge = sim_param['domain_right_edge']
        self.unique_identifier = sim_param['unique_identifier']
        self.cosmological_simulation = sim_param['cosmological_simulation']
        self.current_redshift = sim_param['current_redshift']
        self.omega_lambda = sim_param['omega_lambda']
        self.hubble_constant = sim_param['hubble_constant']
        fileh.close()
        
         
    @classmethod
    def _is_valid(self, *args, **kwargs):
        format = 'Gadget Infrastructure'
        add1 = 'griddded_data_format'
        add2 = 'data_software'
        try:
            fileh = h5py.File(args[0],'r')
            if add1 in fileh['/'].items():
                if add2 in fileh['/'+add1].attrs.keys():
                    if fileh['/'+add1].attrs[add2] == format:
                        fileh.close()
                        return True
            fileh.close()
        except:
            pass
        return False
