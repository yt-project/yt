"""
Python-based grid handler, not to be confused with the SWIG-handler

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
@license:
  Copyright (C) 2007 Matthew Turk.  All Rights Reserved.

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

class EnzoGridBase(EnzoData):
    _spatial = True
    _num_ghost_zones = 0
    """
    Class representing a single Enzo Grid instance
    """
    _grids = None

    def __init__(self, id, filename=None, hierarchy = None):
        """
        Returns an instance of EnzoGrid

        @param hierarchy: EnzoHierarchy, parent hierarchy
        @type hierarchy: L{EnzoHierarchy<EnzoHierarchy>}
        @param id: grid ID (NOT index, which is ID-1)
        @type id: int
        @keyword filename: filename holding grid data
        @type filename: string
        """
        #EnzoData.__init__(self, None, [])
        self.data = {}
        self.field_parameters = {}
        self.fields = []
        self.start_index = None
        self.id = id
        if hierarchy: self.hierarchy = hierarchy
        if filename: self.set_filename(filename)
        self.overlap_masks = [None, None, None]
        self._overlap_grids = [None, None, None]
        self._file_access_pooling = False

    def __len__(self):
        return na.prod(self.ActiveDimensions)

    def _generate_field(self, fieldName):
        """
        Generates, or attempts to generate, a field not found in the data file

        @param fieldName: field to generate
        @type fieldName: string
        """
        if fieldInfo.has_key(fieldName):
            # First we check the validator
            try:
                fieldInfo[fieldName].check_available(self)
            except NeedsGridType, ngt_exception:
                # This is only going to be raised if n_gz > 0
                n_gz = ngt_exception.ghost_zones
                f_gz = ngt_exception.fields
                gz_grid = self.retrieve_ghost_zones(n_gz, f_gz)
                temp_array = fieldInfo[fieldName](gz_grid)
                sl = [slice(n_gz,-n_gz)] * 3
                self[fieldName] = temp_array[sl]
            else:
                self[fieldName] = fieldInfo[fieldName](self)
        else: # Can't find the field, try as it might
            raise exceptions.KeyError, fieldName

    def get_data(self, field):
        """
        Returns a field or set of fields for a key or set of keys
        """
        if not self.data.has_key(field):
            if field in self.hierarchy.field_list:
                conv_factor = 1.0
                if fieldInfo.has_key(field):
                    conv_factor = fieldInfo[field]._convert_function(self)
                self[field] = self.readDataFast(field) * conv_factor
            else:
                self._generate_field(field)
        return self.data[field]

    def clear_all_grid_references(self):
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

    def _prepare_grid(self):
        """
        Copies all the appropriate attributes from the hierarchy
        """
        # This is definitely the slowest part of generating the hierarchy
        # Now we give it pointers to all of its attributes
        # Note that to keep in line with Enzo, we have broken PEP-8
        h = self.hierarchy # cache it
        self.Dimensions = h.gridDimensions[self.id-1]
        self.StartIndices = h.gridStartIndices[self.id-1]
        self.EndIndices = h.gridEndIndices[self.id-1]
        self.LeftEdge = h.gridLeftEdge[self.id-1]
        self.RightEdge = h.gridRightEdge[self.id-1]
        self.Level = h.gridLevels[self.id-1,0]
        self.Time = h.gridTimes[self.id-1,0]
        self.NumberOfParticles = h.gridNumberOfParticles[self.id-1,0]
        self.ActiveDimensions = (self.EndIndices - self.StartIndices + 1)
        self.Children = h.gridTree[self.id-1]
        pID = h.gridReverseTree[self.id-1]
        if pID != None and pID != -1:
            self.Parent = h.grids[pID - 1]
        else:
            self.Parent = None

    def _setup_dx(self):
        # So first we figure out what the index is.  We don't assume
        # that dx=dy=dz , at least here.  We probably do elsewhere.
        self.dx = self.hierarchy.gridDxs[self.id-1,0]
        self.dy = self.hierarchy.gridDys[self.id-1,0]
        self.dz = self.hierarchy.gridDzs[self.id-1,0]
        self.data['dx'] = self.dx
        self.data['dy'] = self.dy
        self.data['dz'] = self.dz
        self._corners = self.hierarchy.gridCorners[:,:,self.id-1]

    def _guess_properties_from_parent(self):
        # Okay, we're going to try to guess
        # We know that our grid boundary occurs on the cell boundary of our
        # parent
        le = self.LeftEdge
        self.dx = self.Parent.dx/2.0
        self.dy = self.Parent.dy/2.0
        self.dz = self.Parent.dz/2.0
        ParentLeftIndex = na.rint((self.LeftEdge-self.Parent.LeftEdge)/self.Parent.dx)
        self.start_index = 2*(ParentLeftIndex + self.Parent.get_global_startindex()).astype('int64')
        self.LeftEdge = self.Parent.LeftEdge + self.Parent.dx * ParentLeftIndex
        self.RightEdge = self.LeftEdge + \
                         self.ActiveDimensions*na.array([self.dx,self.dy,self.dz])
        self.hierarchy.gridDxs[self.id-1,0] = self.dx
        self.hierarchy.gridDys[self.id-1,0] = self.dy
        self.hierarchy.gridDzs[self.id-1,0] = self.dz
        self.hierarchy.gridLeftEdge[self.id-1,:] = self.LeftEdge
        self.hierarchy.gridRightEdge[self.id-1,:] = self.RightEdge
        self.hierarchy.gridCorners[:,:,self.id-1] = na.array([ # Unroll!
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
        self.__child_indices = None
        self._setup_dx()

    def get_global_startindex(self):
        if self.start_index != None:
            return self.start_index
        if self.Parent == None:
            start_index = self.LeftEdge / na.array([self.dx, self.dy, self.dz])
            return na.rint(start_index).astype('int64').ravel()
        pdx = na.array([self.Parent.dx, self.Parent.dy, self.Parent.dz]).ravel()
        start_index = (self.Parent.get_global_startindex()) + \
                       na.rint((self.LeftEdge - self.Parent.LeftEdge)/pdx)
        self.start_index = (start_index*2).astype('int64').ravel()
        return self.start_index

    def _generate_overlap_masks(self, axis, LE, RE):
        """
        Generate a mask that shows which cells overlap with other cells on
        different grids.  (If fed appropriate subsets, can be constrained to
        current level.
        Use algorithm described at http://www.gamedev.net/reference/articles/article735.asp

        @param axis: axis  along which line of sight is drawn
        @type axis: int
        @param LE: LeftEdge positions to check against
        @type LE: array of floats
        @param RE: RightEdge positions to check against
        @type RE: array of floats
        """
        x = x_dict[axis]
        y = y_dict[axis]
        cond1 = self.RightEdge[x] >= LE[:,x]
        cond2 = self.LeftEdge[x] <= RE[:,x]
        cond3 = self.RightEdge[y] >= LE[:,y]
        cond4 = self.LeftEdge[y] <= RE[:,y]
        self.overlap_masks[axis]=na.logical_and(na.logical_and(cond1, cond2), \
                                                na.logical_and(cond3, cond4))
    def __repr__(self):
        return "Grid_%04i" % (self.id)

    def __int__(self):
        return self.id

    def clear_data(self):
        self._del_child_mask()
        self._del_child_indices()
        if hasattr(self, 'coarseData'):
            del self.coarseData
        if hasattr(self, 'retVal'):
            del self.retVal
        EnzoData.clear_data(self)
        self._setup_dx()

    def set_filename(self, filename):
        if self.hierarchy._strip_path:
            self.filename = os.path.join(self.hierarchy.directory,
                                         os.path.basename(filename))
        elif filename[0] == os.path.sep:
            self.filename = filename
        else:
            self.filename = os.path.join(self.hierarchy.directory, filename)
        return

    def find_max(self, field):
        """
        Returns value, coordinate of maximum value in this gird

        @param field: field to check
        @type field: string
        """
        coord=nd.maximum_position(self[field]*self.child_mask)
        val = self[field][coord]
        return val, coord

    def find_min(self, field):
        """
        Returns value, coordinate of minimum value in this gird

        @param field: field to check
        @type field: string
        """
        coord=nd.minimum_position(self[field])
        val = self[field][coord]
        return val, coord

    def get_position(self, coord):
        """
        Returns position of a coordinate

        @param coord: position to check
        @type coord: array of floats
        """
        pos = (coord + 0.0) * self.dx + self.LeftEdge
        # Should 0.0 be 0.5?
        return pos

    def clear_all(self):
        """
        Clears all datafields from memory.
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

    def get_enzo_grid(self):
        """
        This attempts to get an instance of this particular grid from the SWIG
        interface.  Note that it first checks to see if the ParameterFile has
        been instantiated.
        """
        if self.hierarchy.eiTopGrid == None:
            self.hierarchy.initializeEnzoInterface()
        p=re.compile("Grid = %s\n" % (self.id))
        h=open(self.hierarchyFilename,"r").read()
        m=re.search(p,h)
        h=open(self.hierarchyFilename,"r")
        retVal = yt.enki.EnzoInterface.fseek(h, long(m.end()), 0)
        self.eiGrid=yt.enki.EnzoInterface.grid()
        cwd = os.getcwd() # Hate doing this, need to for relative pathnames
        os.chdir(self.hierarchy.directory)
        self.eiGrid.ReadGrid(h, 1)
        os.chdir(cwd)
        mylog.debug("Grid read with SWIG")

    def export_amira(self, filename, fields, timestep = 1, a5Filename=None, gid=0):
        fields = ensure_list(fields)
        deltas = na.array([self.dx,self.dy,self.dz],dtype='float64')
        tn = "time-%i" % (timestep)
        ln = "level-%i" % (self.Level)
        for field in fields:
            iorigin = (self.LeftEdge/deltas).astype('int64')
            new_h5 = tables.openFile(filename % {'field' : field}, "a")
            f = self[field].transpose().reshape(self.ActiveDimensions)
            new_h5.createArray("/","grid-%i" % (self.id), f)
            del f
            node = new_h5.getNode("/","grid-%i" % (self.id))
            node.setAttr("level",self.Level)
            node.setAttr("timestep",timestep)
            node.setAttr("time",self.Time)
            node.setAttr("cctk_bbox",na.array([0,0,0,0,0,0],dtype='int32'))
            node.setAttr("cctk_nghostzones",na.array([0,0,0],dtype='int32'))
            node.setAttr("delta",deltas)
            node.setAttr("origin",self.LeftEdge)
            node.setAttr("iorigin",iorigin*(2**(self.hierarchy.maxLevel - self.Level)))
            new_h5.close()
            if a5Filename != None:
                new_h5 = tables.openFile(a5Filename % {'field' : field}, "a")
                new_h5.createGroup("/%s/%s" % (tn, ln),"grid-%i" % (gid))
                node=new_h5.getNode("/%s/%s" % (tn, ln),"grid-%i" % (gid))
                node._f_setAttr("dims",self.ActiveDimensions)
                node._f_setAttr("ghostzoneFlags",na.array([0,0,0,0,0,0],dtype='int32'))
                node._f_setAttr("integerOrigin",(self.LeftEdge/deltas).astype('int64'))
                node._f_setAttr("numGhostzones",na.array([0,0,0],dtype='int32'))
                node._f_setAttr("origin",self.LeftEdge)
                node._f_setAttr("referenceDataPath","/"+"grid-%i" % (self.id))
                fn = os.path.basename(filename % {'field' : field})
                node._f_setAttr("referenceFileName", fn)
                new_h5.close()

    def _get_projection(self, axis, field, zeroOut, weight=None, func=na.sum):
        """
        Projects along an axis.  Currently in flux.  Shouldn't be called
        directly.
        """
        if weight == None:
            maskedData = self[field].copy()
            weightData = na.ones(maskedData.shape)
        else:
            maskedData = self[field] * self[weight]
            weightData = self[weight].copy()
        if zeroOut:
            maskedData[self.child_indices]=0
            weightData[self.child_indices]=0
            toCombineMask = na.logical_and.reduce(self.child_mask, axis).astype('int64')
        a = {0:self.dx, 1:self.dy, 2:self.dz}
        fullProj = func(maskedData,axis)*a[axis] # Gives correct shape
        weightProj = func(weightData,axis)*a[axis]
        if not zeroOut:
            toCombineMask = na.ones(fullProj.shape, dtype='int64')
        toCombineMask = toCombineMask.astype('int64')
        cmI = na.indices(fullProj.shape)
        xind = cmI[0,:].ravel()
        yind = cmI[1,:].ravel()
        dx = a[x_dict[axis]]
        dy = a[y_dict[axis]]
        xpoints = xind + na.rint(self.LeftEdge[x_dict[axis]]/dx).astype('int64')
        ypoints = yind + na.rint(self.LeftEdge[y_dict[axis]]/dy).astype('int64')
        return [xpoints, ypoints, fullProj.ravel(), toCombineMask.ravel(), weightProj.ravel()]

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

    #@time_execution
    def __generate_child_mask(self):
        """
        Generates self.child_mask, which is zero where child grids exist (and
        thus, where higher resolution data is available.)
        """
        self.__child_mask = na.ones(self.ActiveDimensions, 'int32')
        for child in self.Children:
            # Now let's get our overlap
            si = [None]*3
            ei = [None]*3
            startIndex = na.rint((child.LeftEdge - self.LeftEdge)/self.dx)
            endIndex = na.rint((child.RightEdge - self.LeftEdge)/self.dx)
            for i in xrange(3):
                si[i] = startIndex[i]
                ei[i] = endIndex[i]
            self.__child_mask[si[0]:ei[0], si[1]:ei[1], si[2]:ei[2]] = 0
        self.__child_indices = na.where(self.__child_mask==0)

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
        self['x'], self['y'], self['z'] = (ind+0.5)*self.dx+LE

    __child_mask = None
    __child_indices = None

    child_mask = property(fget=_get_child_mask, fdel=_del_child_mask)
    child_indices = property(fget=_get_child_indices, fdel = _del_child_indices)

    def retrieve_ghost_zones(self, n_zones, fields, all_levels=False):
        # We will attempt this by creating a datacube that is exactly bigger
        # than the grid by nZones*dx in each direction
        new_left_edge = self.LeftEdge - n_zones * self.dx
        new_right_edge = self.RightEdge + n_zones * self.dx
        # Something different needs to be done for the root grid, though
        level = self.Level
        if all_levels:
            level = self.hierarchy.max_level + 1
        cube = self.hierarchy.covering_grid(level,
                        new_left_edge, new_right_edge,
                        self.ActiveDimensions + 2*n_zones, fields,
                        num_ghost_zones=n_zones, use_pbar=False)
        return cube

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
