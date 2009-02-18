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

class AMRGridPatch(AMRData):
    _spatial = True
    _num_ghost_zones = 0
    _grids = None
    _id_offset = 1

    _type_name = 'grid'
    _con_args = ('id', 'filename')

    def __init__(self, id, filename=None, hierarchy = None):
        self.data = {}
        self.field_parameters = {}
        self.fields = []
        self.start_index = None
        self.id = id
        if (id % 1e4) == 0: mylog.debug("Prepared grid %s", id)
        if hierarchy: self.hierarchy = weakref.proxy(hierarchy)
        if filename: self.set_filename(filename)
        self.pf = self.hierarchy.parameter_file # weakref already

    def _generate_field(self, field):
        if self.pf.field_info.has_key(field):
            # First we check the validator
            try:
                self.pf.field_info[field].check_available(self)
            except NeedsGridType, ngt_exception:
                # This is only going to be raised if n_gz > 0
                n_gz = ngt_exception.ghost_zones
                f_gz = ngt_exception.fields
                gz_grid = self.retrieve_ghost_zones(n_gz, f_gz)
                temp_array = self.pf.field_info[field](gz_grid)
                sl = [slice(n_gz,-n_gz)] * 3
                self[field] = temp_array[sl]
            else:
                self[field] = self.pf.field_info[field](self)
        else: # Can't find the field, try as it might
            raise exceptions.KeyError, field
    
    def get_data(self, field):
        """
        Returns a field or set of fields for a key or set of keys
        """
        if not self.data.has_key(field):
            if field in self.hierarchy.field_list:
                conv_factor = 1.0
                if self.pf.field_info.has_key(field):
                    conv_factor = self.pf.field_info[field]._convert_function(self)
                try:
                    if hasattr(self.hierarchy, 'queue'):
                        temp = self.hierarchy.queue.pop(self, field)
                    else:
                        temp = self.readDataFast(field)
                    self[field] = temp * conv_factor
                except self._read_exception, exc:
                    if field in self.pf.field_info:
                        if self.pf.field_info[field].particle_type:
                            self[field] = na.array([],dtype='int64')
                        elif self.pf.field_info[field].not_in_all:
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
        self.dds = na.array([self.hierarchy.gridDxs[id,0],
                                     self.hierarchy.gridDys[id,0],
                                     self.hierarchy.gridDzs[id,0]])
        self.data['dx'], self.data['dy'], self.data['dz'] = self.dds

    @property
    def _corners(self):
        return self.hierarchy.gridCorners[:,:,self.id - self._id_offset]

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
        return "Grid_%04i" % (self.id)

    def __int__(self):
        return self.id

    def clear_all_grid_references(self):
        """
        This clears out all references this grid has to any others, as
        well as the hierarchy.  It's like extra-cleaning after clear_data.
        """
        self.clear_all_derived_quantities()
        if hasattr(self, 'hierarchy'):
            del self.hierarchy
        if hasattr(self, 'Parent'):
            if self.Parent != None:
                self.Parent.clear_all_grid_references()
            del self.Parent
        if hasattr(self, 'Children'):
            for i in self.Children:
                if i != None:
                    del i
            del self.Children

    def clear_data(self):
        """
        Clear out the following things: child_mask, child_indices,
        all fields, all field parameters.
        """
        self._del_child_mask()
        self._del_child_indices()
        if hasattr(self, 'coarseData'):
            del self.coarseData
        if hasattr(self, 'retVal'):
            del self.retVal
        AMRData.clear_data(self)
        self._setup_dx()

    def _prepare_grid(self):
        """
        Copies all the appropriate attributes from the hierarchy
        """
        # This is definitely the slowest part of generating the hierarchy
        # Now we give it pointers to all of its attributes
        # Note that to keep in line with Enzo, we have broken PEP-8
        h = self.hierarchy # cache it
        my_ind = self.id - self._id_offset
        self.ActiveDimensions = h.gridEndIndices[my_ind] \
                              - h.gridStartIndices[my_ind] + 1
        self.LeftEdge = h.gridLeftEdge[my_ind]
        self.RightEdge = h.gridRightEdge[my_ind]
        self.Level = h.gridLevels[my_ind,0]
        # This might be needed for streaming formats
        #self.Time = h.gridTimes[my_ind,0]
        self.NumberOfParticles = h.gridNumberOfParticles[my_ind,0]
        self.Children = h.gridTree[my_ind]
        pID = h.gridReverseTree[my_ind]
        if pID != None and pID != -1:
            self.Parent = weakref.proxy(h.grids[pID - self._id_offset])
        else:
            self.Parent = None

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
        if self.__child_mask != None:
            mylog.warning("Overriding child_mask attribute!  This is probably unwise!")
        self.__child_mask = newCM

    def _set_child_indices(self, newCI):
        if self.__child_indices != None:
            mylog.warning("Overriding child_indices attribute!  This is probably unwise!")
        self.__child_indices = newCI

    def _get_child_mask(self):
        if self.__child_mask == None:
            self.__generate_child_mask()
        return self.__child_mask

    def _get_child_indices(self):
        if self.__child_indices == None:
            self.__generate_child_mask()
        return self.__child_indices

    def _del_child_indices(self):
        try:
            del self.__child_indices
        except AttributeError:
            pass
        self.__child_indices = None

    def _del_child_mask(self):
        try:
            del self.__child_mask
        except AttributeError:
            pass
        self.__child_mask = None

    def _get_child_index_mask(self):
        if self.__child_index_mask is None:
            self.__generate_child_index_mask()
        return self.__child_index_mask

    def _del_child_index_mask(self):
        try:
            del self.__child_index_mask
        except AttributeError:
            pass
        self.__child_index_mask = None

    #@time_execution
    def __fill_child_mask(self, child, mask, tofill):
        startIndex = na.maximum(0, na.rint(
                    (child.LeftEdge - self.LeftEdge)/self.dds))
        endIndex = na.minimum(na.rint(
                    (child.RightEdge - self.LeftEdge)/self.dds),
                              self.ActiveDimensions)
        startIndex = na.maximum(0, startIndex)
        mask[startIndex[0]:endIndex[0],
             startIndex[1]:endIndex[1],
             startIndex[2]:endIndex[2]] = tofill

    def __generate_child_mask(self):
        """
        Generates self.child_mask, which is zero where child grids exist (and
        thus, where higher resolution data is available.)
        """
        self.__child_mask = na.ones(self.ActiveDimensions, 'int32')
        for child in self.Children:
            self.__fill_child_mask(child, self.__child_mask, 0)
        self.__child_indices = (self.__child_mask==0) # bool, possibly redundant

    def __generate_child_index_mask(self):
        """
        Generates self.child_index_mask, which is -1 where there is no child,
        and otherwise has the ID of the grid that resides there.
        """
        self.__child_index_mask = na.zeros(self.ActiveDimensions, 'int32') - 1
        for child in self.Children:
            self.__fill_child_mask(child, self.__child_index_mask,
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

    __child_mask = None
    __child_indices = None
    __child_index_mask = None

    child_mask = property(fget=_get_child_mask, fdel=_del_child_mask)
    child_index_mask = property(fget=_get_child_index_mask, fdel=_del_child_index_mask)
    child_indices = property(fget=_get_child_indices, fdel = _del_child_indices)

    def retrieve_ghost_zones(self, n_zones, fields, all_levels=False,
                             smoothed=False):
        # We will attempt this by creating a datacube that is exactly bigger
        # than the grid by nZones*dx in each direction
        new_left_edge = self.LeftEdge - n_zones * self.dds
        new_right_edge = self.RightEdge + n_zones * self.dds
        # Something different needs to be done for the root grid, though
        level = self.Level
        if all_levels:
            level = self.hierarchy.max_level + 1
        args = (level, new_left_edge, new_right_edge)
        kwargs = {'dims': self.ActiveDimensions + 2*n_zones,
                  'num_ghost_zones':n_zones,
                  'use_pbar':False, 'fields':fields}
        if smoothed:
            cube = self.hierarchy.smoothed_covering_grid(
                level, new_left_edge, new_right_edge, **kwargs)
        else:
            cube = self.hierarchy.covering_grid(
                level, new_left_edge, new_right_edge, **kwargs)
        return cube

    def get_vertex_centered_data(self, field):
        cg = self.retrieve_ghost_zones(1, field, smoothed=True)
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

    def _save_data_state(self):
        self.__current_data_keys = self.data.keys()
        if self.__child_mask != None:
            self.__current_child_mask == True
        else:
            self.__current_child_mask = False

        if self.__child_indices != None:
            self.__current_child_indices == True
        else:
            self.__current_child_indices = False

    def _restore_data_state(self):
        if not self.__current_child_mask:
            self._del_child_mask()
        if not self.__current_child_indices:
            self._del_child_indices()
        for key in data.keys():
            if key not in self.__current_data_keys:
                del self.data[key]

class EnzoGridBase(AMRGridPatch):
    """
    Class representing a single Enzo Grid instance.
    """
    def __init__(self, id, filename=None, hierarchy = None):
        """
        Returns an instance of EnzoGrid with *id*, associated with
        *filename* and *hierarchy*.
        """
        #All of the field parameters will be passed to us as needed.
        AMRGridPatch.__init__(self, id, filename, hierarchy)

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
        self['dx'] = self.Parent['dx']/rf
        self['dy'] = self.Parent['dy']/rf
        self['dz'] = self.Parent['dz']/rf
        ParentLeftIndex = na.rint((self.LeftEdge-self.Parent.LeftEdge)/self.Parent.dds)
        self.start_index = rf*(ParentLeftIndex + self.Parent.get_global_startindex()).astype('int64')
        self.LeftEdge = self.Parent.LeftEdge + self.Parent.dds * ParentLeftIndex
        self.RightEdge = self.LeftEdge + self.ActiveDimensions*self.dds
        self.hierarchy.gridDxs[my_ind,0] = self['dx']
        self.hierarchy.gridDys[my_ind,0] = self['dy']
        self.hierarchy.gridDzs[my_ind,0] = self['dz']
        self.hierarchy.gridLeftEdge[my_ind,:] = self.LeftEdge
        self.hierarchy.gridRightEdge[my_ind,:] = self.RightEdge
        self.hierarchy.gridCorners[:,:,my_ind] = na.array([ # Unroll!
            [self.LeftEdge[0], self.LeftEdge[1], self.LeftEdge[2]],
            [self.RightEdge[0], self.LeftEdge[1], self.LeftEdge[2]],
            [self.RightEdge[0], self.RightEdge[1], self.LeftEdge[2]],
            [self.RightEdge[0], self.RightEdge[1], self.RightEdge[2]],
            [self.LeftEdge[0], self.RightEdge[1], self.RightEdge[2]],
            [self.LeftEdge[0], self.LeftEdge[1], self.RightEdge[2]],
            [self.RightEdge[0], self.LeftEdge[1], self.RightEdge[2]],
            [self.LeftEdge[0], self.RightEdge[1], self.LeftEdge[2]],
            ], dtype='float64')
        self.__child_mask = None
        self.__child_index_mask = None
        self.__child_indices = None
        self._setup_dx()

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

class OrionGridBase(AMRGridPatch):
    _id_offset = 0
    def __init__(self, LeftEdge, RightEdge, index, level, filename, offset, dimensions,start,stop,paranoia=False):
        AMRGridPatch.__init__(self, index)
        self.filename = filename
        self._offset = offset
        self._paranoid = paranoia
        
        # should error check this
        self.ActiveDimensions = dimensions.copy()#.transpose()
        self.start = start.copy()#.transpose()
        self.stop = stop.copy()#.transpose()
        self.LeftEdge  = LeftEdge.copy()
        self.RightEdge = RightEdge.copy()
        self.index = index
        self.Level = level

    def get_global_startindex(self):
        return self.start + na.rint(self.pf["DomainLeftEdge"]/self.dds)

    def _prepare_grid(self):
        """
        Copies all the appropriate attributes from the hierarchy
        """
        # This is definitely the slowest part of generating the hierarchy
        # Now we give it pointers to all of its attributes
        # Note that to keep in line with Enzo, we have broken PEP-8
        h = self.hierarchy # cache it
        self.StartIndices = h.gridStartIndices[self.id]
        self.EndIndices = h.gridEndIndices[self.id]
        h.gridLevels[self.id,0] = self.Level
        h.gridLeftEdge[self.id,:] = self.LeftEdge[:]
        h.gridRightEdge[self.id,:] = self.RightEdge[:]
        self.Time = h.gridTimes[self.id,0]
        self.NumberOfParticles = h.gridNumberOfParticles[self.id,0]
        self.Children = h.gridTree[self.id]
        pIDs = h.gridReverseTree[self.id]
        if len(pIDs) > 0:
            self.Parent = [weakref.proxy(h.grids[pID]) for pID in pIDs]
        else:
            self.Parent = []

