"""
Data structures for Gadget.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: UCSD
Author: Chris Moody <cemoody@ucsc.edu>
Affiliation: UCSC
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2007-2010 Matthew Turk.  All Rights Reserved.

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
import numpy as np
import ipdb

from yt.funcs import *
from yt.data_objects.grid_patch import \
    AMRGridPatch
from yt.data_objects.hierarchy import \
    AMRHierarchy
from yt.data_objects.static_output import \
    StaticOutput

from .fields import GadgetFieldContainer

class GadgetGrid(AMRGridPatch):
    _id_offset = 0
    def __init__(self, hierarchy, id,dimensions,left_index,
            level,parent_id,particle_count,Parent=[],Children=[]):
        AMRGridPatch.__init__(self, id, filename = hierarchy.filename,
                              hierarchy = hierarchy)
        self.id = id
        self.address = '/data/grid_%010i'%id
        self.dimensions = dimensions
        self.left_index = left_index
        self.level = level
        self.Level = level
        self.parent_id = parent_id
        self.particle_count = particle_count
        
        try:
            padd =self.address+'/particles'
            self.particle_types = self.hierarchy._handle[padd].keys()
        except:
            self.particle_types =  ()
        self.filename = hierarchy.filename
        
        self.Parent = Parent
        self.Children = Children
        
    def _setup_dx(self):
        # So first we figure out what the index is.  We don't assume
        # that dx=dy=dz , at least here.  We probably do elsewhere.
        if len(self.Parent)>0:
            self.dds = self.Parent[0].dds / self.pf.refine_by
        else:
            LE, RE = self.hierarchy.grid_left_edge[self.id,:], \
                     self.hierarchy.grid_right_edge[self.id,:]
            self.dds = np.array((RE-LE)/self.ActiveDimensions)
        if self.pf.dimensionality < 2: self.dds[1] = 1.0
        if self.pf.dimensionality < 3: self.dds[2] = 1.0
        self.data['dx'], self.data['dy'], self.data['dz'] = self.dds
        
    def __repr__(self):
        return 'GadgetGrid_%05i'%self.id
        
class GadgetHierarchy(AMRHierarchy):
    grid = GadgetGrid

    def __init__(self, pf, data_style='gadget_hdf5'):
        self.field_info = GadgetFieldContainer()
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
        npa = np.array
        
        grid_levels = npa(f['/grid_level'],dtype='int')
        self.grid_left_edge  = npa(f['/grid_left_index'],dtype='int')
        self.grid_dimensions = npa(f['/grid_dimensions'],dtype='int')
        self.grid_right_edge = self.grid_left_edge + self.grid_dimensions 
        grid_particle_count = npa(f['/grid_particle_count'],dtype='int')
        self.grid_parent_id = npa(f['/grid_parent_id'],dtype='int')
        self.max_level = np.max(self.grid_levels)
        
        args = zip(range(self.num_grids),grid_levels,self.grid_parent_id,
            self.grid_left_edge,self.grid_dimensions, grid_particle_count)
        grids = [self.grid(self,j,d,le,lvl,p,n) for j,lvl,p, le, d, n in args]
        self.grids = npa(grids,dtype='object')
        
        
        
    def _populate_grid_objects(self):    
        ipdb.set_trace()
        for g in self.grids:
            if g.parent_id>=0:
                parent = self.grids[g.parent_id]
                g.Parent = g.Parent + [parent,]
                parent.Children = parent.Children + [g,]
        for g in self.grids:
            g._prepare_grid()
            g._setup_dx()
            
        
    def _setup_unknown_fields(self):
        pass

    def _setup_derived_fields(self):
        self.derived_field_list = []

class GadgetStaticOutput(StaticOutput):
    _hierarchy_class = GadgetHierarchy
    _fieldinfo_class = GadgetFieldContainer
    def __init__(self, filename,storage_filename=None) :
        self.storage_filename = storage_filename
        self.field_info = self._fieldinfo_class()
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
        fileh = h5py.File(args[0],'r')
        if add1 in fileh['/'].items():
            if add2 in fileh['/'+add1].attrs.keys():
                if fileh['/'+add1].attrs[add2] == format:
                    return True
        return False
