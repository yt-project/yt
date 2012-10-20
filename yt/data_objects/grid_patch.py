"""
Python-based grid handler, not to be confused with the SWIG-handler

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
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

import exceptions
import pdb
import weakref

import numpy as np

from yt.funcs import *
from yt.utilities.definitions import x_dict, y_dict

from yt.data_objects.data_containers import YTFieldData
from yt.utilities.definitions import x_dict, y_dict
from .field_info_container import \
    NeedsGridType, \
    NeedsOriginalGrid, \
    NeedsDataField, \
    NeedsProperty, \
    NeedsParameter
from yt.geometry.selection_routines import convert_mask_to_indices

class AMRGridPatch(object):
    _spatial = True
    _num_ghost_zones = 0
    _grids = None
    _id_offset = 1

    _type_name = 'grid'
    _skip_add = True
    _con_args = ('id', 'filename')
    OverlappingSiblings = None

    __slots__ = ['field_data', 'field_parameters', 'id', 'hierarchy', 'pf',
                 'ActiveDimensions', 'LeftEdge', 'RightEdge', 'Level',
                 'NumberOfParticles', 'Children', 'Parent',
                 'start_index', 'filename', '__weakref__', 'dds',
                 '_child_mask', '_child_indices', '_child_index_mask',
                 '_parent_id', '_children_ids',
                 '_last_mask', '_last_selector_id']

    def __init__(self, id, filename=None, hierarchy=None):
        self.field_data = YTFieldData()
        self.field_parameters = {}
        self.id = id
        if hierarchy: self.hierarchy = weakref.proxy(hierarchy)
        self.pf = self.hierarchy.parameter_file  # weakref already
        self._child_mask = self._child_indices = self._child_index_mask = None
        self.start_index = None
        self._last_mask = None
        self._last_selector_id = None

    def get_global_startindex(self):
        """
        Return the integer starting index for each dimension at the current
        level.

        """
        if self.start_index is not None:
            return self.start_index
        if self.Parent == None:
            left = self.LeftEdge - self.pf.domain_left_edge
            start_index = left / self.dds
            return np.rint(start_index).astype('int64').ravel()

        pdx = self.Parent.dds
        start_index = (self.Parent.get_global_startindex()) + \
                       np.rint((self.LeftEdge - self.Parent.LeftEdge) / pdx)
        self.start_index = (start_index * self.pf.refine_by).astype('int64').ravel()
        return self.start_index

    def get_field_parameter(self, name, default=None):
        """
        This is typically only used by derived field functions, but it returns
        parameters used to generate fields.

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
        """ Checks if a field parameter is set. """
        return self.field_parameters.has_key(name)

    def convert(self, datatype):
        """
        This will attempt to convert a given unit to cgs from code units. It
        either returns the multiplicative factor or throws a KeyError.

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
                if f_gz is None:
                    f_gz = self.pf.field_info[field].get_dependencies(
                            pf = self.pf).requested
                gz_grid = self.retrieve_ghost_zones(n_gz, f_gz, smoothed=True)
                temp_array = self.pf.field_info[field](gz_grid)
                sl = [slice(n_gz, -n_gz)] * 3
                self[field] = temp_array[sl]
            else:
                self[field] = self.pf.field_info[field](self)
        else: # Can't find the field, try as it might
            raise exceptions.KeyError(field)

    def has_key(self, key):
        return (key in self.field_data)

    def __getitem__(self, key):
        """
        Returns a single field.  Will add if necessary.
        """
        if not self.field_data.has_key(key):
            self.get_data(key)
        return self.field_data[key]

    def __setitem__(self, key, val):
        """
        Sets a field to be some other value.
        """
        self.field_data[key] = val

    def __delitem__(self, key):
        """
        Deletes a field
        """
        del self.field_data[key]

    def keys(self):
        return self.field_data.keys()

    def get_data(self, field, convert = True):
        """
        Returns a field or set of fields for a key or set of keys
        """
        if not self.field_data.has_key(field):
            if field in self.hierarchy.field_list:
                conv_factor = 1.0
                if self.pf.field_info.has_key(field) and convert == True:
                    conv_factor = self.pf.field_info[field]._convert_function(self)
                if self.pf.field_info[field].particle_type and \
                   self.NumberOfParticles == 0:
                    # because this gets upcast to float
                    self[field] = np.array([],dtype='int64')
                    return self.field_data[field]
                try:
                    temp = self.hierarchy.io.pop(self, field)
                    self[field] = np.multiply(temp, conv_factor, temp)
                except self.hierarchy.io._read_exception, exc:
                    if field in self.pf.field_info:
                        if self.pf.field_info[field].not_in_all:
                            self[field] = np.zeros(self.ActiveDimensions, dtype='float64')
                        else:
                            raise
                    else: raise
            else:
                self._generate_field(field)
        return self.field_data[field]

    def _setup_dx(self):
        # So first we figure out what the index is.  We don't assume
        # that dx=dy=dz, at least here.  We probably do elsewhere.
        id = self.id - self._id_offset
        if self.Parent is not None:
            self.dds = self.Parent.dds / self.pf.refine_by
        else:
            LE, RE = self.hierarchy.grid_left_edge[id,:], \
                     self.hierarchy.grid_right_edge[id,:]
            self.dds = np.array((RE - LE) / self.ActiveDimensions)
        if self.pf.dimensionality < 2: self.dds[1] = self.pf.domain_right_edge[1] - self.pf.domain_left_edge[1]
        if self.pf.dimensionality < 3: self.dds[2] = self.pf.domain_right_edge[2] - self.pf.domain_left_edge[2]
        self.field_data['dx'], self.field_data['dy'], self.field_data['dz'] = self.dds

    @property
    def _corners(self):
        return np.array([ # Unroll!
            [self.LeftEdge[0],  self.LeftEdge[1],  self.LeftEdge[2]],
            [self.RightEdge[0], self.LeftEdge[1],  self.LeftEdge[2]],
            [self.RightEdge[0], self.RightEdge[1], self.LeftEdge[2]],
            [self.RightEdge[0], self.RightEdge[1], self.RightEdge[2]],
            [self.LeftEdge[0],  self.RightEdge[1], self.RightEdge[2]],
            [self.LeftEdge[0],  self.LeftEdge[1],  self.RightEdge[2]],
            [self.RightEdge[0], self.LeftEdge[1],  self.RightEdge[2]],
            [self.LeftEdge[0],  self.RightEdge[1], self.LeftEdge[2]],
            ], dtype='float64')

    @property
    def shape(self):
        return tuple(self.ActiveDimensions)

    def _generate_overlap_masks(self, axis, LE, RE):
        """
        Generate a mask that shows which cells overlap with arbitrary arrays
        *LE* and *RE*) of edges, typically grids, along *axis*.
        Use algorithm described at http://www.gamedev.net/reference/articles/article735.asp

        """
        x = x_dict[axis]
        y = y_dict[axis]
        cond = self.RightEdge[x] >= LE[:,x]
        cond = np.logical_and(cond, self.LeftEdge[x] <= RE[:,x])
        cond = np.logical_and(cond, self.RightEdge[y] >= LE[:,y])
        cond = np.logical_and(cond, self.LeftEdge[y] <= RE[:,y])
        return cond

    def __repr__(self):
        return "AMRGridPatch_%04i" % (self.id)

    def __int__(self):
        return self.id

    def clear_data(self):
        """
        Clear out the following things: child_mask, child_indices, all fields,
        all field parameters.

        """
        self._del_child_mask()
        self._del_child_indices()
        self.field_data.clear()
        self._setup_dx()

    def check_child_masks(self):
        return self._child_mask, self._child_indices

    def _prepare_grid(self):
        """ Copies all the appropriate attributes from the hierarchy. """
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
        self.NumberOfParticles = h.grid_particle_count[my_ind, 0]

    def __len__(self):
        return np.prod(self.ActiveDimensions)

    def find_max(self, field):
        """ Returns value, index of maximum value of *field* in this grid. """
        coord1d = (self[field] * self.child_mask).argmax()
        coord = np.unravel_index(coord1d, self[field].shape)
        val = self[field][coord]
        return val, coord

    def find_min(self, field):
        """ Returns value, index of minimum value of *field* in this grid. """
        coord1d = (self[field] * self.child_mask).argmin()
        coord = np.unravel_index(coord1d, self[field].shape)
        val = self[field][coord]
        return val, coord

    def get_position(self, index):
        """ Returns center position of an *index*. """
        pos = (index + 0.5) * self.dds + self.LeftEdge
        return pos

    def clear_all(self):
        """
        Clears all datafields from memory and calls
        :meth:`clear_derived_quantities`.

        """
        for key in self.keys():
            del self.field_data[key]
        del self.field_data
        if hasattr(self,"retVal"):
            del self.retVal
        self.field_data = YTFieldData()
        self.clear_derived_quantities()

    def clear_derived_quantities(self):
        """ Clears coordinates, child_indices, child_mask. """
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
    def __fill_child_mask(self, child, mask, tofill, dlevel = 1):
        rf = self.pf.refine_by
        if dlevel != 1:
            rf = rf**dlevel
        gi, cgi = self.get_global_startindex(), child.get_global_startindex()
        startIndex = np.maximum(0, cgi / rf - gi)
        endIndex = np.minimum((cgi + child.ActiveDimensions) / rf - gi,
                              self.ActiveDimensions)
        endIndex += (startIndex == endIndex)
        mask[startIndex[0]:endIndex[0],
             startIndex[1]:endIndex[1],
             startIndex[2]:endIndex[2]] = tofill

    def __generate_child_mask(self):
        """
        Generates self.child_mask, which is zero where child grids exist (and
        thus, where higher resolution data is available).

        """
        self._child_mask = np.ones(self.ActiveDimensions, 'bool')
        for child in self.Children:
            self.__fill_child_mask(child, self._child_mask, 0)
        if self.OverlappingSiblings is not None:
            for sibling in self.OverlappingSiblings:
                self.__fill_child_mask(sibling, self._child_mask, 0)
        
        self._child_indices = (self._child_mask==0) # bool, possibly redundant

    def __generate_child_index_mask(self):
        """
        Generates self.child_index_mask, which is -1 where there is no child,
        and otherwise has the ID of the grid that resides there.

        """
        self._child_index_mask = np.zeros(self.ActiveDimensions, 'int32') - 1
        for child in self.Children:
            self.__fill_child_mask(child, self._child_index_mask,
                                   child.id)
        if self.OverlappingSiblings is not None:
            for sibling in self.OverlappingSiblings:
                self.__fill_child_mask(sibling, self._child_index_mask,
                                       sibling.id)

    def _get_coords(self):
        if self.__coords == None: self._generate_coords()
        return self.__coords

    def _set_coords(self, new_c):
        if self.__coords != None:
            mylog.warning("Overriding coords attribute!  This is probably unwise!")
        self.__coords = new_c

    def _del_coords(self):
        del self.__coords
        self.__coords = None

    def _generate_coords(self):
        """
        Creates self.coords, which is of dimensions (3, ActiveDimensions)

        """
        ind = np.indices(self.ActiveDimensions)
        left_shaped = np.reshape(self.LeftEdge, (3, 1, 1, 1))
        self['x'], self['y'], self['z'] = (ind + 0.5) * self.dds + left_shaped

    child_mask = property(fget=_get_child_mask, fdel=_del_child_mask)
    child_index_mask = property(fget=_get_child_index_mask, fdel=_del_child_index_mask)
    child_indices = property(fget=_get_child_indices, fdel = _del_child_indices)

    def retrieve_ghost_zones(self, n_zones, fields, all_levels=False,
                             smoothed=False):
        # We will attempt this by creating a datacube that is exactly bigger
        # than the grid by nZones*dx in each direction
        nl = self.get_global_startindex() - n_zones
        nr = nl + self.ActiveDimensions + 2 * n_zones
        new_left_edge = nl * self.dds + self.pf.domain_left_edge
        new_right_edge = nr * self.dds + self.pf.domain_left_edge

        # Something different needs to be done for the root grid, though
        level = self.Level
        if all_levels:
            level = self.hierarchy.max_level + 1
        args = (level, new_left_edge, new_right_edge)
        kwargs = {'dims': self.ActiveDimensions + 2*n_zones,
                  'num_ghost_zones':n_zones,
                  'use_pbar':False, 'fields':fields}
        # This should update the arguments to set the field parameters to be
        # those of this grid.
        kwargs.update(self.field_parameters)
        if smoothed:
            cube = self.hierarchy.smoothed_covering_grid(
                level, new_left_edge, **kwargs)
        else:
            cube = self.hierarchy.covering_grid(level, new_left_edge, **kwargs)

        return cube

    def get_vertex_centered_data(self, field, smoothed=True, no_ghost=False):
        new_field = np.zeros(self.ActiveDimensions + 1, dtype='float64')

        if no_ghost:
            of = self[field]
            new_field[:-1,:-1,:-1] += of
            new_field[:-1,:-1,1:] += of
            new_field[:-1,1:,:-1] += of
            new_field[:-1,1:,1:] += of
            new_field[1:,:-1,:-1] += of
            new_field[1:,:-1,1:] += of
            new_field[1:,1:,:-1] += of
            new_field[1:,1:,1:] += of
            np.multiply(new_field, 0.125, new_field)
            if self.pf.field_info[field].take_log:
                new_field = np.log10(new_field)

            new_field[:,:, -1] = 2.0*new_field[:,:,-2] - new_field[:,:,-3]
            new_field[:,:, 0]  = 2.0*new_field[:,:,1] - new_field[:,:,2]
            new_field[:,-1, :] = 2.0*new_field[:,-2,:] - new_field[:,-3,:]
            new_field[:,0, :]  = 2.0*new_field[:,1,:] - new_field[:,2,:]
            new_field[-1,:,:] = 2.0*new_field[-2,:,:] - new_field[-3,:,:]
            new_field[0,:,:]  = 2.0*new_field[1,:,:] - new_field[2,:,:]

            if self.pf.field_info[field].take_log:
                np.power(10.0, new_field, new_field)
        else:
            cg = self.retrieve_ghost_zones(1, field, smoothed=smoothed)
            np.add(new_field, cg[field][1: ,1: ,1: ], new_field)
            np.add(new_field, cg[field][:-1,1: ,1: ], new_field)
            np.add(new_field, cg[field][1: ,:-1,1: ], new_field)
            np.add(new_field, cg[field][1: ,1: ,:-1], new_field)
            np.add(new_field, cg[field][:-1,1: ,:-1], new_field)
            np.add(new_field, cg[field][1: ,:-1,:-1], new_field)
            np.add(new_field, cg[field][:-1,:-1,1: ], new_field)
            np.add(new_field, cg[field][:-1,:-1,:-1], new_field)
            np.multiply(new_field, 0.125, new_field)

        return new_field

    def icoords(self, dobj):
        mask = self.select(dobj.selector)
        if mask is None: return np.empty((0,3), dtype='int64')
        coords = convert_mask_to_indices(mask, mask.sum())
        coords += self.get_global_startindex()[None, :]
        return coords

    def fcoords(self, dobj):
        mask = self.select(dobj.selector)
        if mask is None: return np.empty((0,3), dtype='float64')
        coords = convert_mask_to_indices(mask, mask.sum()).astype("float64")
        coords += 0.5
        coords *= self.dds[None, :]
        coords += self.LeftEdge[None, :]
        return coords

    def fwidth(self, dobj):
        mask = self.select(dobj.selector)
        if mask is None: return np.empty((0,3), dtype='float64')
        coords = np.empty((mask.sum(), 3), dtype='float64')
        for axis in range(3):
            coords[:,axis] = self.dds[axis]
        return coords

    def ires(self, dobj):
        mask = self.select(dobj.selector)
        if mask is None: return np.empty(0, dtype='int64')
        coords = np.empty(mask.sum(), dtype='int64')
        coords[:] = self.Level
        return coords

    def select(self, selector):
        if id(selector) == self._last_selector_id:
            return self._last_mask
        self._last_mask = selector.fill_mask(self)
        self._last_selector_id = id(selector)
        return self._last_mask

    def count(self, selector):
        if id(selector) == self._last_selector_id:
            if self._last_mask is None: return 0
            return self._last_mask.sum()
        self.select(selector)
        return self.count(selector)

    def count_particles(self, selector, x, y, z):
        # We don't cache the selector results
        count = selector.count_points(x,y,z)
        return count

    def select_particles(self, selector, x, y, z):
        mask = selector.select_points(x,y,z)
        return mask
