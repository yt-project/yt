"""
Python-based grid handler, not to be confused with the SWIG-handler

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2007-2009 Matthew Turk.  All Rights Reserved.

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

from yt.lagos import *
#import yt.enki, gc
from yt.funcs import *

class AMRGridPatch(object):
    _spatial = True
    _num_ghost_zones = 0
    _grids = None
    _id_offset = 1

    _type_name = 'grid'
    _skip_add = True
    _con_args = ('id', 'filename')

    __slots__ = ['data', 'field_parameters', 'id', 'hierarchy', 'pf',
                 'ActiveDimensions', 'LeftEdge', 'RightEdge', 'Level',
                 'NumberOfParticles', 'Children', 'Parent',
                 'start_index', 'filename', '__weakref__', 'dds',
                 '_child_mask', '_child_indices', '_child_index_mask',
                 '_parent_id', '_children_ids']
    def __init__(self, id, filename = None, hierarchy = None):
        self.data = {}
        self.field_parameters = {}
        self.id = id
        if hierarchy: self.hierarchy = weakref.proxy(hierarchy)
        self.pf = self.hierarchy.parameter_file # weakref already
        self._child_mask = self._child_indices = self._child_index_mask = None
        self.start_index = None

    def get_global_startindex(self):
        """
        Return the integer starting index for each dimension at the current
        level.
        """
        if self.start_index != None:
            return self.start_index
        if self.Parent == None:
            start_index = self.LeftEdge / self.dds
            return na.rint(start_index).astype('int64').ravel()
        pdx = self.Parent.dds
        start_index = (self.Parent.get_global_startindex()) + \
                       na.rint((self.LeftEdge - self.Parent.LeftEdge)/pdx)
        self.start_index = (start_index*self.pf["RefineBy"]).astype('int64').ravel()
        return self.start_index


    def get_field_parameter(self, name, default=None):
        """
        This is typically only used by derived field functions, but
        it returns parameters used to generate fields.
        """
        if self.field_parameters.has_key(name):
            return self.field_parameters[name]
        else:
            return default

    def set_field_parameter(self, name, val):
        """
        Here we set up dictionaries that get passed up and down and ultimately
        to derived fields.
        """
        self.field_parameters[name] = val

    def has_field_parameter(self, name):
        """
        Checks if a field parameter is set.
        """
        return self.field_parameters.has_key(name)

    def convert(self, datatype):
        """
        This will attempt to convert a given unit to cgs from code units.
        It either returns the multiplicative factor or throws a KeyError.
        """
        return self.pf[datatype]

    def __repr__(self):
        # We'll do this the slow way to be clear what's going on
        s = "%s (%s): " % (self.__class__.__name__, self.pf)
        s += ", ".join(["%s=%s" % (i, getattr(self,i))
                       for i in self._con_args])
        return s

    def _generate_field(self, field):
        if self.pf.field_info.has_key(field):
            # First we check the validator
            try:
                self.pf.field_info[field].check_available(self)
            except NeedsGridType, ngt_exception:
                # This is only going to be raised if n_gz > 0
                n_gz = ngt_exception.ghost_zones
                f_gz = ngt_exception.fields
                gz_grid = self.retrieve_ghost_zones(n_gz, f_gz, smoothed=True)
                temp_array = self.pf.field_info[field](gz_grid)
                sl = [slice(n_gz,-n_gz)] * 3
                self[field] = temp_array[sl]
            else:
                self[field] = self.pf.field_info[field](self)
        else: # Can't find the field, try as it might
            raise exceptions.KeyError, field

    def has_key(self, key):
        return (key in self.data)

    def __getitem__(self, key):
        """
        Returns a single field.  Will add if necessary.
        """
        if not self.data.has_key(key):
            self.get_data(key)
        return self.data[key]

    def __setitem__(self, key, val):
        """
        Sets a field to be some other value.
        """
        self.data[key] = val

    def __delitem__(self, key):
        """
        Deletes a field
        """
        del self.data[key]

    def keys(self):
        return self.data.keys()
    
    def get_data(self, field):
        """
        Returns a field or set of fields for a key or set of keys
        """
        if not self.data.has_key(field):
            if field in self.hierarchy.field_list:
                conv_factor = 1.0
                if self.pf.field_info.has_key(field):
                    conv_factor = self.pf.field_info[field]._convert_function(self)
                if self.pf.field_info[field].particle_type and \
                   self.NumberOfParticles == 0:
                    # because this gets upcast to float
                    self[field] = na.array([],dtype='int64')
                    return self.data[field]
                try:
                    temp = self.hierarchy.io.pop(self, field)
                    self[field] = na.multiply(temp, conv_factor, temp)
                except self.hierarchy.io._read_exception, exc:
                    if field in self.pf.field_info:
                        if self.pf.field_info[field].not_in_all:
                            self[field] = na.zeros(self.ActiveDimensions, dtype='float64')
                        else:
                            raise
                    else: raise
            else:
                self._generate_field(field)
        return self.data[field]

    def _setup_dx(self):
        # So first we figure out what the index is.  We don't assume
        # that dx=dy=dz , at least here.  We probably do elsewhere.
        id = self.id - self._id_offset
        if self.Parent is not None:
            self.dds = self.Parent.dds / self.pf["RefineBy"]
        else:
            LE, RE = self.hierarchy.grid_left_edge[id,:], \
                     self.hierarchy.grid_right_edge[id,:]
            self.dds = na.array((RE-LE)/self.ActiveDimensions)
        if self.pf["TopGridRank"] < 2: self.dds[1] = 1.0
        if self.pf["TopGridRank"] < 3: self.dds[2] = 1.0
        self.data['dx'], self.data['dy'], self.data['dz'] = self.dds

    @property
    def _corners(self):
        return na.array([ # Unroll!
            [self.LeftEdge[0],  self.LeftEdge[1],  self.LeftEdge[2]],
            [self.RightEdge[0], self.LeftEdge[1],  self.LeftEdge[2]],
            [self.RightEdge[0], self.RightEdge[1], self.LeftEdge[2]],
            [self.RightEdge[0], self.RightEdge[1], self.RightEdge[2]],
            [self.LeftEdge[0],  self.RightEdge[1], self.RightEdge[2]],
            [self.LeftEdge[0],  self.LeftEdge[1],  self.RightEdge[2]],
            [self.RightEdge[0], self.LeftEdge[1],  self.RightEdge[2]],
            [self.LeftEdge[0],  self.RightEdge[1], self.LeftEdge[2]],
            ], dtype='float64')

    def _generate_overlap_masks(self, axis, LE, RE):
        """
        Generate a mask that shows which cells overlap with arbitrary arrays
        *LE* and *RE*) of edges, typically grids, along *axis*.
        Use algorithm described at http://www.gamedev.net/reference/articles/article735.asp
        """
        x = x_dict[axis]
        y = y_dict[axis]
        cond = self.RightEdge[x] >= LE[:,x]
        cond = na.logical_and(cond, self.LeftEdge[x] <= RE[:,x])
        cond = na.logical_and(cond, self.RightEdge[y] >= LE[:,y])
        cond = na.logical_and(cond, self.LeftEdge[y] <= RE[:,y])
        return cond
   
    def __repr__(self):
        return "AMRGridPatch_%04i" % (self.id)

    def __int__(self):
        return self.id

    def clear_data(self):
        """
        Clear out the following things: child_mask, child_indices,
        all fields, all field parameters.
        """
        self._del_child_mask()
        self._del_child_indices()
        self.data.clear()
        self._setup_dx()

    def check_child_masks(self):
        return self._child_mask, self._child_indices

    def _prepare_grid(self):
        """
        Copies all the appropriate attributes from the hierarchy
        """
        # This is definitely the slowest part of generating the hierarchy
        # Now we give it pointers to all of its attributes
        # Note that to keep in line with Enzo, we have broken PEP-8
        h = self.hierarchy # cache it
        my_ind = self.id - self._id_offset
        self.ActiveDimensions = h.grid_dimensions[my_ind]
        self.LeftEdge = h.grid_left_edge[my_ind]
        self.RightEdge = h.grid_right_edge[my_ind]
        h.grid_levels[my_ind, 0] = self.Level
        # This might be needed for streaming formats
        #self.Time = h.gridTimes[my_ind,0]
        self.NumberOfParticles = h.grid_particle_count[my_ind,0]

    def __len__(self):
        return na.prod(self.ActiveDimensions)

    def find_max(self, field):
        """
        Returns value, index of maximum value of *field* in this gird
        """
        coord1d=(self[field]*self.child_mask).argmax()
        coord=na.unravel_index(coord1d, self[field].shape)
        val = self[field][coord]
        return val, coord

    def find_min(self, field):
        """
        Returns value, index of minimum value of *field* in this gird
        """
        coord1d=(self[field]*self.child_mask).argmin()
        coord=na.unravel_index(coord1d, self[field].shape)
        val = self[field][coord]
        return val, coord

    def get_position(self, index):
        """
        Returns center position of an *index*
        """
        pos = (index + 0.5) * self.dds + self.LeftEdge
        return pos

    def clear_all(self):
        """
        Clears all datafields from memory and calls
        :meth:`clear_derived_quantities`.
        """
        for key in self.keys():
            del self.data[key]
        del self.data
        if hasattr(self,"retVal"):
            del self.retVal
        self.data = {}
        self.clear_derived_quantities()

    def clear_derived_quantities(self):
        """
        Clears coordinates, child_indices, child_mask.
        """
        # Access the property raw-values here
        del self.child_mask
        del self.child_ind

    def _set_child_mask(self, newCM):
        if self._child_mask != None:
            mylog.warning("Overriding child_mask attribute!  This is probably unwise!")
        self._child_mask = newCM

    def _set_child_indices(self, newCI):
        if self._child_indices != None:
            mylog.warning("Overriding child_indices attribute!  This is probably unwise!")
        self._child_indices = newCI

    def _get_child_mask(self):
        if self._child_mask == None:
            self.__generate_child_mask()
        return self._child_mask

    def _get_child_indices(self):
        if self._child_indices == None:
            self.__generate_child_mask()
        return self._child_indices

    def _del_child_indices(self):
        try:
            del self._child_indices
        except AttributeError:
            pass
        self._child_indices = None

    def _del_child_mask(self):
        try:
            del self._child_mask
        except AttributeError:
            pass
        self._child_mask = None

    def _get_child_index_mask(self):
        if self._child_index_mask is None:
            self.__generate_child_index_mask()
        return self._child_index_mask

    def _del_child_index_mask(self):
        try:
            del self._child_index_mask
        except AttributeError:
            pass
        self._child_index_mask = None

    #@time_execution
    def __fill_child_mask(self, child, mask, tofill):
        rf = self.pf["RefineBy"]
        gi, cgi = self.get_global_startindex(), child.get_global_startindex()
        startIndex = na.maximum(0, cgi/rf - gi)
        endIndex = na.minimum( (cgi+child.ActiveDimensions)/rf - gi,
                              self.ActiveDimensions)
        endIndex += (startIndex == endIndex)
        mask[startIndex[0]:endIndex[0],
             startIndex[1]:endIndex[1],
             startIndex[2]:endIndex[2]] = tofill

    def __generate_child_mask(self):
        """
        Generates self.child_mask, which is zero where child grids exist (and
        thus, where higher resolution data is available.)
        """
        self._child_mask = na.ones(self.ActiveDimensions, 'int32')
        for child in self.Children:
            self.__fill_child_mask(child, self._child_mask, 0)
        self._child_indices = (self._child_mask==0) # bool, possibly redundant

    def __generate_child_index_mask(self):
        """
        Generates self.child_index_mask, which is -1 where there is no child,
        and otherwise has the ID of the grid that resides there.
        """
        self._child_index_mask = na.zeros(self.ActiveDimensions, 'int32') - 1
        for child in self.Children:
            self.__fill_child_mask(child, self._child_index_mask,
                                   child.id)

    def _get_coords(self):
        if self.__coords == None: self._generate_coords()
        return self.__coords

    def _set_coords(self, newC):
        if self.__coords != None:
            mylog.warning("Overriding coords attribute!  This is probably unwise!")
        self.__coords = newC

    def _del_coords(self):
        del self.__coords
        self.__coords = None

    def _generate_coords(self):
        """
        Creates self.coords, which is of dimensions (3,ActiveDimensions)
        """
        #print "Generating coords"
        ind = na.indices(self.ActiveDimensions)
        LE = na.reshape(self.LeftEdge,(3,1,1,1))
        self['x'], self['y'], self['z'] = (ind+0.5)*self.dds+LE

    child_mask = property(fget=_get_child_mask, fdel=_del_child_mask)
    child_index_mask = property(fget=_get_child_index_mask, fdel=_del_child_index_mask)
    child_indices = property(fget=_get_child_indices, fdel = _del_child_indices)

    def retrieve_ghost_zones(self, n_zones, fields, all_levels=False,
                             smoothed=False):
        # We will attempt this by creating a datacube that is exactly bigger
        # than the grid by nZones*dx in each direction
        nl = self.get_global_startindex() - n_zones
        nr = nl + self.ActiveDimensions + 2*n_zones
        new_left_edge = nl * self.dds + self.pf["DomainLeftEdge"]
        new_right_edge = nr * self.dds + self.pf["DomainLeftEdge"]
        # Something different needs to be done for the root grid, though
        level = self.Level
        if all_levels:
            level = self.hierarchy.max_level + 1
        args = (level, new_left_edge, new_right_edge)
        kwargs = {'dims': self.ActiveDimensions + 2*n_zones,
                  'num_ghost_zones':n_zones,
                  'use_pbar':False, 'fields':fields}
        if smoothed:
            #cube = self.hierarchy.smoothed_covering_grid(
            #    level, new_left_edge, new_right_edge, **kwargs)
            cube = self.hierarchy.si_covering_grid(
                level, new_left_edge, **kwargs)
        else:
            cube = self.hierarchy.covering_grid(
                level, new_left_edge, **kwargs)
        return cube

    def get_vertex_centered_data(self, field, smoothed=True):
        cg = self.retrieve_ghost_zones(1, field, smoothed=smoothed)
        # Bounds should be cell-centered
        bds = na.array(zip(cg.left_edge+cg.dds/2.0, cg.right_edge-cg.dds/2.0)).ravel()
        interp = TrilinearFieldInterpolator(na.log10(cg[field]), bds, ['x','y','z'])
        ad = self.ActiveDimensions + 1
        x,y,z = na.mgrid[self.LeftEdge[0]:self.RightEdge[0]:ad[0]*1j,
                         self.LeftEdge[1]:self.RightEdge[1]:ad[1]*1j,
                         self.LeftEdge[2]:self.RightEdge[2]:ad[2]*1j]
        dd = {'x':x,'y':y,'z':z}
        scalars = 10**interp(dict(x=x,y=y,z=z))
        return scalars

class EnzoGrid(AMRGridPatch):
    """
    Class representing a single Enzo Grid instance.
    """

    __slots__ = []
    def __init__(self, id, hierarchy):
        """
        Returns an instance of EnzoGrid with *id*, associated with
        *filename* and *hierarchy*.
        """
        #All of the field parameters will be passed to us as needed.
        AMRGridPatch.__init__(self, id, filename = None, hierarchy = hierarchy)
        self._children_ids = []
        self._parent_id = -1
        self.Level = -1

    def _guess_properties_from_parent(self):
        """
        We know that our grid boundary occurs on the cell boundary of our
        parent.  This can be a very expensive process, but it is necessary
        in some hierarchys, where yt is unable to generate a completely
        space-filling tiling of grids, possibly due to the finite accuracy in a
        standard Enzo hierarchy file.
        """
        rf = self.pf["RefineBy"]
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
        """
        Intelligently set the filename.
        """
        if self.hierarchy._strip_path:
            self.filename = os.path.join(self.hierarchy.directory,
                                         os.path.basename(filename))
        elif filename[0] == os.path.sep:
            self.filename = filename
        else:
            self.filename = os.path.join(self.hierarchy.directory, filename)
        return

    def __repr__(self):
        return "EnzoGrid_%04i" % (self.id)

    @property
    def Parent(self):
        if self._parent_id == -1: return None
        return self.hierarchy.grids[self._parent_id - self._id_offset]

    @property
    def Children(self):
        return [self.hierarchy.grids[cid - self._id_offset]
                for cid in self._children_ids]

class EnzoGridInMemory(EnzoGrid):
    __slots__ = ['proc_num']
    def set_filename(self, filename):
        pass

class OrionGrid(AMRGridPatch):
    _id_offset = 0
    def __init__(self, LeftEdge, RightEdge, index, level, filename, offset, dimensions,start,stop,paranoia=False,**kwargs):
        AMRGridPatch.__init__(self, index,**kwargs)
        self.filename = filename
        self._offset = offset
        self._paranoid = paranoia
        
        # should error check this
        self.ActiveDimensions = (dimensions.copy()).astype('int32')#.transpose()
        self.start_index = start.copy()#.transpose()
        self.stop_index = stop.copy()#.transpose()
        self.LeftEdge  = LeftEdge.copy()
        self.RightEdge = RightEdge.copy()
        self.index = index
        self.Level = level

    def get_global_startindex(self):
        return self.start_index

    def _prepare_grid(self):
        """
        Copies all the appropriate attributes from the hierarchy
        """
        # This is definitely the slowest part of generating the hierarchy
        # Now we give it pointers to all of its attributes
        # Note that to keep in line with Enzo, we have broken PEP-8
        h = self.hierarchy # cache it
        #self.StartIndices = h.gridStartIndices[self.id]
        #self.EndIndices = h.gridEndIndices[self.id]
        h.grid_levels[self.id,0] = self.Level
        h.grid_left_edge[self.id,:] = self.LeftEdge[:]
        h.grid_right_edge[self.id,:] = self.RightEdge[:]
        #self.Time = h.gridTimes[self.id,0]
        #self.NumberOfParticles = h.gridNumberOfParticles[self.id,0]
        self.field_indexes = h.field_indexes
        self.Children = h.gridTree[self.id]
        pIDs = h.gridReverseTree[self.id]
        if len(pIDs) > 0:
            self.Parent = [weakref.proxy(h.grids[pID]) for pID in pIDs]
        else:
            self.Parent = None

    def _setup_dx(self):
        # So first we figure out what the index is.  We don't assume
        # that dx=dy=dz , at least here.  We probably do elsewhere.
        id = self.id - self._id_offset
        if self.Parent is not None:
            self.dds = self.Parent[0].dds / self.pf["RefineBy"]
        else:
            LE, RE = self.hierarchy.grid_left_edge[id,:], \
                     self.hierarchy.grid_right_edge[id,:]
            self.dds = na.array((RE-LE)/self.ActiveDimensions)
        if self.pf["TopGridRank"] < 2: self.dds[1] = 1.0
        if self.pf["TopGridRank"] < 3: self.dds[2] = 1.0
        self.data['dx'], self.data['dy'], self.data['dz'] = self.dds

    def __repr__(self):
        return "OrionGrid_%04i" % (self.id)

class ProtoGadgetGrid(object):
    def __init__(self, level, left_edge, right_edge, particles, parent = None):
        # generate the left/right edge from the octant and the level argument
        self.level = level
        self.left_edge = left_edge
        self.right_edge = right_edge
        self.particles = particles
        self.parent = parent # We break convention here -- this is a protogrid
        self.children = []

    def split(self):
        #split_axis = (self.left_edge + self.right_edge) / 2.0
        W = 0.5 * (self.right_edge - self.left_edge) # This is the half-width
        keep = na.ones(self.particles.shape[1], dtype='bool')
        for i in xrange(2):
            for j in xrange(2):
                for k in xrange(2):
                    LE = self.left_edge.copy()
                    LE += na.array([i, j, k], dtype='float64') * W
                    RE = LE + W
                    pi = self.select_particles(LE, RE)
                    keep &= (~pi)
                    g = ProtoGadgetGrid(self.level + 1, LE, RE,
                                        self.particles[:,pi], self)
                    self.children.append(g)
                    yield g
        self.particles = self.particles[:,keep]

    def select_particles(self, LE, RE):
        pi = na.all( (self.particles >= LE[:,None])
                   & (self.particles <  RE[:,None]), axis=0)
        return pi

    def refine(self, npart = 128):
        gs = [self]
        if self.particles.shape[1] < npart: return gs
        for grid in self.split():
            gs += grid.refine(npart)
        return gs

class GadgetGrid(AMRGridPatch):

    _id_offset = 0

    def __init__(self, id, hierarchy, proto_grid):
        AMRGridPatch.__init__(self, id, hierarchy = hierarchy)
        self.Children = []
        if proto_grid.parent is None:
            self.Parent = None
        else:
            self.Parent = proto_grid.parent.real_grid
            self.Parent.Children.append(self)
        self['particle_position_x'] = proto_grid.particles[0,:]
        self['particle_position_y'] = proto_grid.particles[1,:]
        self['particle_position_z'] = proto_grid.particles[2,:]
        self.NumberOfParticles = proto_grid.particles.shape[1]
        self.Level = proto_grid.level
        self.LeftEdge = proto_grid.left_edge
        self.RightEdge = proto_grid.right_edge
        # Our dx is a bit fluid here, so we defer
        dims = self.pf["TopGridDimensions"]
        # Hard code to refineby 2
        dds = 1.0 / (dims * 2**self.Level)
        ad = na.rint((self.RightEdge - self.LeftEdge) / dds)
        self.ActiveDimensions = ad.astype('int64')

    def __repr__(self):
        return "GadgetGrid_%04i" % (self.id)

