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

    def __init__(self, id, filename, hierarchy, **kwargs):
        AMRGridPatch.__init__(self, id, hierarchy = hierarchy)
        self.id = id
        self.filename = filename
        self.Children = [] #grid objects
        self.Parent = None
        self.Level = 0
        self.LeftEdge = [0,0,0]
        self.RightEdge = [0,0,0]
        self.IsLeaf = False
        self.N = 0
        self.Address = ''
        self.NumberOfParticles = self.N
        self.ActiveDimensions = na.array([0,0,0])
        self._id_offset = 0
        self.start_index = na.array([0,0,0])
        
        for key,val in kwargs.items():
            if key in dir(self):
                #if it's one of the predefined values
                setattr(self,key,val)
        
    #def __repr__(self):
    #    return "GadgetGrid_%04i" % (self.Address)
    
    def get_global_startindex(self):
        return self.start_index
        
    def _prepare_grid(self):
        #all of this info is already included in the snapshots
        pass
        #h = self.hierarchy
        #h.grid_levels[self.Address,0]=self.Level
        #h.grid_left_edge[self.Address,:]=self.LeftEdge[:]
        #h.grid_right_edge[self.Address,:]=self.RightEdge[:]
    
    def _setup_dx(self):
        # So first we figure out what the index is.  We don't assume
        # that dx=dy=dz , at least here.  We probably do elsewhere.
        id = self.id
        LE, RE = self.LeftEdge,self.RightEdge
        self.dds = na.array((RE-LE)/self.ActiveDimensions)
        if self.pf.dimensionality < 2: self.dds[1] = 1.0
        if self.pf.dimensionality < 3: self.dds[2] = 1.0
        self.data['dx'], self.data['dy'], self.data['dz'] = self.dds

class GadgetHierarchy(AMRHierarchy):
    grid = GadgetGrid

    def __init__(self, pf, data_style='gadget_hdf5'):
        self.field_info = GadgetFieldContainer()
        self.directory = os.path.dirname(pf.parameter_filename)
        self.data_style = data_style
        self._handle = h5py.File(pf.parameter_filename)
        AMRHierarchy.__init__(self, pf, data_style)
        self._handle.close()

    def _initialize_data_storage(self):
        pass

    def _detect_fields(self):
        #example string:
        #"(S'VEL'\np1\nS'ID'\np2\nS'MASS'\np3\ntp4\n."
        #fields are surrounded with '
        fields_string=self._handle['root'].attrs['fieldnames']
        #splits=fields_string.split("'")
        #pick out the odd fields
        #fields= [splits[j] for j in range(1,len(splits),2)]
        self.field_list = cPickle.loads(fields_string)
    
    def _setup_classes(self):
        dd = self._get_data_reader_dict()
        AMRHierarchy._setup_classes(self, dd)
        self.object_types.sort()

    def _count_grids(self):
        fh = self._handle #shortcut
        #nodes in the hdf5 file are the same as grids
        #in yt
        #the num of levels and total nodes is already saved
        self._levels   = self.pf._get_param('maxlevel')
        self.num_grids = self.pf._get_param('numnodes')
        
    def _parse_hierarchy(self):
        #for every box, define a self.grid(level,edge1,edge2) 
        #with particle counts, dimensions
        f = self._handle #shortcut
        
        root = f['root']
        grids,numnodes = self._walk_nodes(None,root,[])
        dims = [self.pf.max_grid_size for grid in grids]
        LE = [grid.LeftEdge for grid in grids]
        RE = [grid.RightEdge for grid in grids]
        levels = [grid.Level for grid in grids]
        counts = [(grid.N if grid.IsLeaf else 0) for grid in grids]
        self.grids = na.array(grids,dtype='object')
        self.grid_dimensions[:] = na.array(dims, dtype='int64')
        self.grid_left_edge[:] = na.array(LE, dtype='float64')
        self.grid_right_edge[:] = na.array(RE, dtype='float64')
        self.grid_levels.flat[:] = na.array(levels, dtype='int32')
        self.grid_particle_count.flat[:] = na.array(counts, dtype='int32')
            
    def _walk_nodes(self,parent,node,grids,idx=0):
        pi = cPickle.loads
        loc = node.attrs['h5address']
        
        kwargs = {}
        kwargs['Address'] = loc
        kwargs['Parent'] = parent
        kwargs['Axis']  = self.pf._get_param('divideaxis',location=loc)
        kwargs['Level']  = self.pf._get_param('level',location=loc)
        kwargs['LeftEdge'] = self.pf._get_param('leftedge',location=loc) 
        kwargs['RightEdge'] = self.pf._get_param('rightedge',location=loc)
        kwargs['IsLeaf'] = self.pf._get_param('isleaf',location=loc)
        kwargs['N'] = self.pf._get_param('n',location=loc)
        kwargs['NumberOfParticles'] = self.pf._get_param('n',location=loc)
        dx = self.pf._get_param('dx',location=loc)
        dy = self.pf._get_param('dy',location=loc)
        dz = self.pf._get_param('dz',location=loc)
        divdims = na.array([1,1,1])
        if not kwargs['IsLeaf']: 
            divdims[kwargs['Axis']] = 2
        kwargs['ActiveDimensions'] = divdims
        #Active dimensions:
        #This is the number of childnodes, along with dimensiolaity
        #ie, binary tree can be (2,1,1) but octree is (2,2,2)
        
        idx+=1
        #pdb.set_trace()
        children = []
        if not kwargs['IsLeaf']:
            for child in node.values():
                children,idx=self._walk_nodes(node,child,children,idx=idx)
        
        kwargs['Children'] = children
        grid = self.grid(idx,self.pf.parameter_filename,self,**kwargs)
        grids += children
        grids += [grid,]
        return grids,idx

    def _populate_grid_objects(self):
        for g in self.grids:
            g._prepare_grid()
        self.max_level = cPickle.loads(self._handle['root'].attrs['maxlevel'])
    
    def _setup_unknown_fields(self):
        pass

    def _setup_derived_fields(self):
        self.derived_field_list = []

    def _get_grid_children(self, grid):
        #given a grid, use it's address to find subchildren
        pass

class GadgetHierarchyOld(AMRHierarchy):
    #Kept here to compare for the time being
    grid = GadgetGrid

    def __init__(self, pf, data_style):
        self.directory = pf.fullpath
        self.data_style = data_style
        AMRHierarchy.__init__(self, pf, data_style)

    def _count_grids(self):
        # We actually construct our octree here!
        # ...but we do read in our particles, it seems.
        LE = na.zeros(3, dtype='float64')
        RE = na.ones(3, dtype='float64')
        base_grid = ProtoGadgetGrid(0, LE, RE, self.pf.particles)
        self.proto_grids = base_grid.refine(8)
        self.num_grids = len(self.proto_grids)
        self.max_level = max( (g.level for g in self.proto_grids) )

    def _setup_classes(self):
        dd = self._get_data_reader_dict()
        AMRHierarchy._setup_classes(self, dd)
        self.object_types.sort()

    def _parse_hierarchy(self):
        grids = []
        # We need to fill in dims, LE, RE, level, count
        dims, LE, RE, levels, counts = [], [], [], [], []
        self.proto_grids.sort(key = lambda a: a.level)
        for i, pg in enumerate(self.proto_grids):
            g = self.grid(i, self, pg)
            pg.real_grid = g
            grids.append(g)
            dims.append(g.ActiveDimensions)
            LE.append(g.LeftEdge)
            RE.append(g.RightEdge)
            levels.append(g.Level)
            counts.append(g.NumberOfParticles)
        del self.proto_grids
        self.grids = na.array(grids, dtype='object')
        self.grid_dimensions[:] = na.array(dims, dtype='int64')
        self.grid_left_edge[:] = na.array(LE, dtype='float64')
        self.grid_right_edge[:] = na.array(RE, dtype='float64')
        self.grid_levels.flat[:] = na.array(levels, dtype='int32')
        self.grid_particle_count.flat[:] = na.array(counts, dtype='int32')

    def _populate_grid_objects(self):
        # We don't need to do anything here
        for g in self.grids: g._setup_dx()

    def _detect_fields(self):
        self.field_list = ['particle_position_%s' % ax for ax in 'xyz']

    def _setup_unknown_fields(self):
        pass

    def _setup_derived_fields(self):
        self.derived_field_list = []

class GadgetStaticOutput(StaticOutput):
    _hierarchy_class = GadgetHierarchy
    _fieldinfo_class = GadgetFieldContainer
    def __init__(self, h5filename,storage_filename=None) :
        StaticOutput.__init__(self, h5filename, 'gadget_hdf5')
        self.storage_filename = storage_filename #Don't know what this is
        self.field_info = self._fieldinfo_class()
        x = self._get_param('maxlevel')**2
        self.max_grid_size = (x,x,x)
        self.parameters["InitialTime"] = 0.0
        # These should be explicitly obtained from the file, but for now that
        # will wait until a reorganization of the source tree and better
        # generalization.
        self.parameters["TopGridRank"] = 3
        self.parameters["RefineBy"] = 2
        self.parameters["DomainLeftEdge"] = self.leftedge
        self.parameters["DomainRightEdge"] = self.rightedge
        
        
    def _parse_parameter_file(self):
        # read the units in from the hdf5 file 
        #fill in self.units dict
        #fill in self.time_units dict (keys: 'days','years', '1')
        
        #import all of the parameter file params 
        #this is NOT originally from the gadget snapshot but instead
        #from the paramfile starting the sim
        skips = ('TITLE','CLASS','VERSION') #these are just hdf5 crap
        fh = h5py.File(self.parameter_filename)
        for kw in fh['root'].attrs.keys():
            if any([skip in kw for skip in skips]):
                continue
            val = fh['root'].attrs[kw]
            if type(val)==type(''):
                try:    val = cPickle.loads(val)
                except: pass
            #also, includes unit info
            setattr(self,kw,val)
            
    def _get_param(self,kw,location='/root'):
        fh = h5py.File(self.parameter_filename)
        val = fh[location].attrs[kw]
        try:    val = cPickle.loads(val)
        except: pass
        return val
            
    def _set_units(self):
        #check out the unit params from _parse_parameter_file and use them
        #code below is all filler
        self.units = {}
        self.time_units = {}
        self.conversion_factors = defaultdict(lambda: 1.0)
        self.time_units['1'] = 1
        self.units['1'] = 1.0
        self.units['unitary'] = 1.0
        self.units['cm'] = 1.0
        seconds = 1 #self["Time"]
        self.time_units['years'] = seconds / (365*3600*24.0)
        self.time_units['days']  = seconds / (3600*24.0)
        for key in yt2orionFieldsDict:
            self.conversion_factors[key] = 1.0
        
        
    @classmethod
    def _is_valid(cls, *args, **kwargs):
        # check for a /root to exist in the h5 file
        try:
            h5f=h5py.File(self.h5filename)
            valid = 'root' in h5f.items()[0]
            h5f.close()
            return valid
        except:
            pass
        return False

