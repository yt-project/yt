"""
Various non-grid data containers.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Author: Britton Smith <Britton.Smith@colorado.edu>
Affiliation: University of Colorado at Boulder
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

def restore_grid_state(func):
    """
    A decorator that takes a function with the API of (self, grid, field)
    and ensures that after the function is called, the field_parameters will
    be returned to normal.
    """
    def save_state(self, grid, field=None):
        old_params = grid.field_parameters
        grid.field_parameters = self.field_parameters
        tr = func(self, grid, field)
        grid.field_parameters = old_params
        return tr
    return save_state

def cache_mask(func):
    """
    For computationally intensive indexing operations, we can cache
    between calls.
    """
    def check_cache(self, grid):
        if isinstance(grid, FakeGridForParticles):
            return func(self, grid)
        elif grid.id not in self._cut_masks:
            cm = func(self, grid)
            self._cut_masks[grid.id] = cm
        return self._cut_masks[grid.id]
    return check_cache

def cache_point_indices(func):
    """
    For computationally intensive indexing operations, we can cache
    between calls.
    """
    def check_cache(self, grid, use_child_mask=True):
        if isinstance(grid, FakeGridForParticles):
            return func(self, grid, use_child_mask)
        elif grid.id not in self._point_indices:
            cm = func(self, grid, use_child_mask)
            self._point_indices[grid.id] = cm
        return self._point_indices[grid.id]
    return check_cache

class FakeGridForParticles(object):
    """
    Mock up a grid to insert particle positions and radii
    into for purposes of confinement in an :class:`AMR3DData`.
    """
    def __init__(self, grid):
        self._corners = grid._corners
        self.field_parameters = {}
        self.data = {'x':grid['particle_position_x'],
                     'y':grid['particle_position_y'],
                     'z':grid['particle_position_z'],
                     'dx':grid['dx'],
                     'dy':grid['dy'],
                     'dz':grid['dz']}
        self.dds = grid.dds.copy()
        self.real_grid = grid
        self.child_mask = 1
        self.ActiveDimensions = self.data['x'].shape
    def __getitem__(self, field):
        if field not in self.data.keys():
            if field == "RadiusCode":
                center = self.field_parameters['center']
                tr = na.sqrt( (self['x'] - center[0])**2.0 +
                              (self['y'] - center[1])**2.0 +
                              (self['z'] - center[2])**2.0 )
            else:
                raise KeyError(field)
        else: tr = self.data[field]
        return tr

class AMRData:
    """
    Generic AMRData container.  By itself, will attempt to
    generate field, read fields (method defined by derived classes)
    and deal with passing back and forth field parameters.
    """
    _grids = None
    _num_ghost_zones = 0
    _con_args = ()

    def __init__(self, pf, fields, **kwargs):
        """
        Typically this is never called directly, but only due to inheritance.
        It associates a :class:`~yt.lagos.StaticOutput` with the class,
        sets its initial set of fields, and the remainder of the arguments
        are passed as field_parameters.
        """
        if pf != None:
            self.pf = pf
            self.hierarchy = pf.hierarchy
        if fields == None: fields = []
        self.fields = ensure_list(fields)[:]
        self.data = {}
        self.field_parameters = {}
        self.__set_default_field_parameters()
        self._cut_masks = {}
        self._point_indices = {}
        for key, val in kwargs.items():
            self.set_field_parameter(key, val)

    def __set_default_field_parameters(self):
        self.set_field_parameter("center",na.zeros(3,dtype='float64'))
        self.set_field_parameter("bulk_velocity",na.zeros(3,dtype='float64'))

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
        return self.hierarchy[datatype]

    def clear_data(self):
        """
        Clears out all data from the AMRData instance, freeing memory.
        """
        for key in self.data.keys():
            del self.data[key]
        del self.data
        if self._grids is not None:
            for grid in self._grids: grid.clear_data()
        self.data = {}

    def has_key(self, key):
        """
        Checks if a data field already exists.
        """
        return self.data.has_key(key)

    def _refresh_data(self):
        """
        Wipes data and rereads/regenerates it from the self.fields.
        """
        self.clear_data()
        self.get_data()

    def __getitem__(self, key):
        """
        Returns a single field.  Will add if necessary.
        """
        if not self.data.has_key(key):
            if key not in self.fields:
                self.fields.append(key)
            self.get_data(key)
        return self.data[key]

    def __setitem__(self, key, val):
        """
        Sets a field to be some other value.
        """
        if key not in self.fields: self.fields.append(key)
        self.data[key] = val

    def __delitem__(self, key):
        """
        Deletes a field
        """
        try:
            del self.fields[self.fields.index(key)]
        except ValueError:
            pass
        del self.data[key]

    def _generate_field_in_grids(self, fieldName):
        pass

    _key_fields = None
    def write_out(self, filename, fields=None, format="%0.16e"):
        if fields is None: fields=sorted(self.data.keys())
        if self._key_fields is None: raise ValueError
        field_order = self._key_fields[:]
        for field in field_order: self[field]
        field_order += [field for field in fields if field not in field_order]
        fid = open(filename,"w")
        fid.write("\t".join(["#"] + field_order + ["\n"]))
        field_data = na.array([self.data[field] for field in field_order])
        for line in range(field_data.shape[1]):
            field_data[:,line].tofile(fid, sep="\t", format=format)
            fid.write("\n")
        fid.close()

    def save_object(self, name, filename = None):
        if filename is not None:
            ds = shelve.open(filename, protocol=-1)
            if name in ds:
                mylog.info("Overwriting %s in %s", name, filename)
            ds[name] = self
            ds.close()
        else:
            self.hierarchy.save_object(self, name)

    def __reduce__(self):
        args = tuple([self.pf._hash(), self._type_name] +
                     [getattr(self, n) for n in self._con_args] +
                     [self.field_parameters])
        return (_reconstruct_object, args)

    def __repr__(self):
        # We'll do this the slow way to be clear what's going on
        s = "%s (%s): " % (self.__class__.__name__, self.pf)
        s += ", ".join(["%s=%s" % (i, getattr(self,i))
                       for i in self._con_args])
        return s

class GridPropertiesMixin(object):

    def select_grids(self, level):
        """
        Return all grids on a given level.
        """
        grids = [g for g in self._grids if g.Level == level]
        return grids

    def __get_levelIndices(self):
        if self.__levelIndices: return self.__levelIndices
        # Otherwise, generate
        # We only have to do this once, so it's not terribly expensive:
        ll = {}
        for level in range(MAXLEVEL):
            t = [i for i in range(len(self._grids)) if self._grids[i].Level == level]
            ll[level] = na.array(t)
        self.__levelIndices = ll
        return self.__levelIndices

    def __set_levelIndices(self, val):
        self.__levelIndices = val

    def __del_levelIndices(self):
        del self.__levelIndices
        self.__levelIndices = None

    __levelIndices = None
    levelIndices = property(__get_levelIndices, __set_levelIndices,
                            __del_levelIndices)

    def __get_gridLeftEdge(self):
        if self.__gridLeftEdge == None:
            self.__gridLeftEdge = na.array([g.LeftEdge for g in self._grids])
        return self.__gridLeftEdge

    def __del_gridLeftEdge(self):
        del self.__gridLeftEdge
        self.__gridLeftEdge = None

    def __set_gridLeftEdge(self, val):
        self.__gridLeftEdge = val

    __gridLeftEdge = None
    gridLeftEdge = property(__get_gridLeftEdge, __set_gridLeftEdge,
                              __del_gridLeftEdge)

    def __get_gridRightEdge(self):
        if self.__gridRightEdge == None:
            self.__gridRightEdge = na.array([g.RightEdge for g in self._grids])
        return self.__gridRightEdge

    def __del_gridRightEdge(self):
        del self.__gridRightEdge
        self.__gridRightEdge = None

    def __set_gridRightEdge(self, val):
        self.__gridRightEdge = val

    __gridRightEdge = None
    gridRightEdge = property(__get_gridRightEdge, __set_gridRightEdge,
                             __del_gridRightEdge)

    def __get_gridLevels(self):
        if self.__gridLevels == None:
            self.__gridLevels = na.array([g.Level for g in self._grids])
        return self.__gridLevels

    def __del_gridLevels(self):
        del self.__gridLevels
        self.__gridLevels = None

    def __set_gridLevels(self, val):
        self.__gridLevels = val

    __gridLevels = None
    gridLevels = property(__get_gridLevels, __set_gridLevels,
                             __del_gridLevels)


class AMR1DData(AMRData, GridPropertiesMixin):
    _spatial = False
    def __init__(self, pf, fields, **kwargs):
        AMRData.__init__(self, pf, fields, **kwargs)
        self._grids = None

    def _generate_field_in_grids(self, field, num_ghost_zones=0):
        for grid in self._grids:
            temp = grid[field]

    def _generate_field(self, field):
        if self.pf.field_info.has_key(field):
            # First we check the validator
            try:
                self.pf.field_info[field].check_available(self)
            except NeedsGridType, ngt_exception:
                # We leave this to be implementation-specific
                self._generate_field_in_grids(field, ngt_exception.ghost_zones)
                return False
            else:
                self[field] = self.pf.field_info[field](self)
                return True
        else: # Can't find the field, try as it might
            raise exceptions.KeyError(field)

    def get_data(self, fields=None, in_grids=False):
        if self._grids == None:
            self._get_list_of_grids()
        points = []
        #if not self.has_key('dx'):
            #self._generate_coords()
        if not fields:
            fields_to_get = self.fields
        else:
            fields_to_get = ensure_list(fields)
        mylog.debug("Going to obtain %s", fields_to_get)
        for field in fields_to_get:
            if self.data.has_key(field):
                continue
            mylog.info("Getting field %s from %s", field, len(self._grids))
            if field not in self.hierarchy.field_list and not in_grids:
                if self._generate_field(field):
                    continue # True means we already assigned it
            self[field] = na.concatenate(
                [self._get_data_from_grid(grid, field)
                 for grid in self._grids])

class AMROrthoRayBase(AMR1DData):
    _key_fields = ['x','y','z','dx','dy','dz']
    _type_name = "ortho_ray"
    _con_args = ('axis', 'coords')
    def __init__(self, axis, coords, fields=None, pf=None, **kwargs):
        """
        Dimensionality is reduced to one, and an ordered list of points at an
        (x,y) tuple along *axis* are available.
        """
        AMR1DData.__init__(self, pf, fields, **kwargs)
        self.axis = axis
        self.px_ax = x_dict[self.axis]
        self.py_ax = y_dict[self.axis]
        self.px_dx = 'd%s'%(axis_names[self.px_ax])
        self.py_dx = 'd%s'%(axis_names[self.py_ax])
        self.px, self.py = coords
        self._refresh_data()

    def _get_list_of_grids(self):
        # This bugs me, but we will give the tie to the LeftEdge
        y = na.where( (self.px >=  self.pf.hierarchy.gridLeftEdge[:,self.px_ax])
                    & (self.px < self.pf.hierarchy.gridRightEdge[:,self.px_ax])
                    & (self.py >=  self.pf.hierarchy.gridLeftEdge[:,self.py_ax])
                    & (self.py < self.pf.hierarchy.gridRightEdge[:,self.py_ax]))
        self._grids = self.hierarchy.grids[y]

    def _get_data_from_grid(self, grid, field):
        # We are orthogonal, so we can feel free to make assumptions
        # for the sake of speed.
        if grid.id not in self._cut_masks:
            gdx = just_one(grid[self.px_dx])
            gdy = just_one(grid[self.py_dx])
            x_coord = int((self.px - grid.LeftEdge[self.px_ax])/gdx)
            y_coord = int((self.py - grid.LeftEdge[self.py_ax])/gdy)
            sl = [None,None,None]
            sl[self.px_ax] = slice(x_coord,x_coord+1,None)
            sl[self.py_ax] = slice(y_coord,y_coord+1,None)
            sl[self.axis] = slice(None)
            self._cut_masks[grid.id] = sl
        else:
            sl = self._cut_masks[grid.id]
        if not iterable(grid[field]):
            gf = grid[field] * na.ones(grid.child_mask[sl].shape)
        else:
            gf = grid[field][sl]
        return gf[na.where(grid.child_mask[sl])]

class AMRRayBase(AMR1DData):
    _type_name = "ray"
    _con_args = ('start_point', 'end_point')
    def __init__(self, start_point, end_point, fields=None, pf=None, **kwargs):
        """
        We accept a start point and an end point and then get all the data
        between those two.
        """
        mylog.warning("This code is poorly tested.  It may give bad data!")
        AMR1DData.__init__(self, pf, fields, **kwargs)
        self.start_point = na.array(start_point)
        self.end_point = na.array(end_point)
        self.vec = self.end_point - self.start_point
        self.center = self.start_point
        self.set_field_parameter('center', self.start_point)
        #self._refresh_data()

    def _get_list_of_grids(self):
        # Get the value of the line at each LeftEdge and RightEdge
        LE = self.pf.h.gridLeftEdge
        RE = self.pf.h.gridRightEdge
        p = na.zeros(self.pf.h.num_grids, dtype='bool')
        # Check left faces first
        for i in range(3):
            i1 = (i+1) % 3
            i2 = (i+2) % 3
            vs = self._get_line_at_coord(LE[:,i], i)
            p = p | ( ( (LE[:,i1] < vs[:,i1]) & (RE[:,i1] > vs[:,i1]) ) \
                    & ( (LE[:,i2] < vs[:,i2]) & (RE[:,i2] > vs[:,i2]) ) )
            vs = self._get_line_at_coord(RE[:,i], i)
            p = p | ( ( (LE[:,i1] < vs[:,i1]) & (RE[:,i1] > vs[:,i1]) ) \
                    & ( (LE[:,i2] < vs[:,i2]) & (RE[:,i2] > vs[:,i2]) ) )
        p = p | ( na.all( LE < self.start_point, axis=1 ) 
                & na.all( RE > self.start_point, axis=1 ) )
        p = p | ( na.all( LE < self.end_point,   axis=1 ) 
                & na.all( RE > self.end_point,   axis=1 ) )
        self._grids = self.hierarchy.grids.copy()#[p]

    def _get_line_at_coord(self, v, index):
        # t*self.vec + self.start_point = self.end_point
        t = (v - self.start_point[index])/self.vec[index]
        t = t.reshape((t.shape[0],1))
        return self.start_point + t*self.vec

    def _get_data_from_grid(self, grid, field):
        mask = na.logical_and(self._get_cut_mask(grid),
                              grid.child_mask)
        return grid[field][mask]
        
    @cache_mask
    def _get_cut_mask(self, grid):
        mask = na.zeros(grid.ActiveDimensions, dtype='int')
        import RTIntegrator as RT
        RT.VoxelTraversal(mask, grid.LeftEdge, grid.RightEdge,
                          grid.dds, self.center, self.vec)
        return mask

class AMR2DData(AMRData, GridPropertiesMixin, ParallelAnalysisInterface):
    _key_fields = ['px','py','pdx','pdy']
    """
    Class to represent a set of :class:`AMRData` that's 2-D in nature, and
    thus does not have as many actions as the 3-D data types.
    """
    _spatial = False
    def __init__(self, axis, fields, pf=None, **kwargs):
        """
        Prepares the AMR2DData, normal to *axis*.  If *axis* is 4, we are not
        aligned with any axis.
        """
        self.axis = axis
        AMRData.__init__(self, pf, fields, **kwargs)
        self.set_field_parameter("axis",axis)
        
    def _convert_field_name(self, field):
        return field

    #@time_execution
    def get_data(self, fields = None):
        """
        Iterates over the list of fields and generates/reads them all.
        """
        # We get it for the values in fields and coords
        # We take a 3-tuple of the coordinate we want to slice through, as well
        # as the axis we're slicing along
        self._get_list_of_grids()
        if not self.has_key('pdx'):
            self._generate_coords()
        if fields == None:
            fields_to_get = self.fields
        else:
            fields_to_get = ensure_list(fields)
        temp_data = {}
        for field in fields_to_get:
            if self.data.has_key(field): continue
            if field not in self.hierarchy.field_list:
                if self._generate_field(field):
                    continue # A "True" return means we did it
            # To ensure that we use data from this object as much as possible,
            # we're going to have to set the same thing several times
            temp_data[field] = na.concatenate(
                [self._get_data_from_grid(grid, field)
                 for grid in self._get_grids()])
            # Now the next field can use this field
            self[field] = temp_data[field] 
        # We finalize
        temp_data = self._mpi_catdict(temp_data)
        # And set, for the next group
        for field in temp_data.keys():
            self[field] = temp_data[field]


    def _generate_field(self, field):
        if self.pf.field_info.has_key(field):
            # First we check the validator
            try:
                self.pf.field_info[field].check_available(self)
            except NeedsGridType, ngt_exception:
                # We leave this to be implementation-specific
                self._generate_field_in_grids(field, ngt_exception.ghost_zones)
                return False
            else:
                self[field] = self.pf.field_info[field](self)
                return True
        else: # Can't find the field, try as it might
            raise exceptions.KeyError(field)

    def _generate_field_in_grids(self, field, num_ghost_zones=0):
        for grid in self._grids:
            temp = grid[field]

    def interpolate_discretize(self, LE, RE, field, side, log_spacing=True):
        """
        This returns a uniform grid of points between *LE* and *RE*,
        interpolated using the nearest neighbor method, with *side* points on a
        side.
        """
        import yt.raven.delaunay as de
        if log_spacing:
            zz = na.log10(self[field])
        else:
            zz = self[field]
        xi, yi = na.array( \
                 na.mgrid[LE[0]:RE[0]:side*1j, \
                          LE[1]:RE[1]:side*1j], 'float64')
        zi = de.Triangulation(self['px'],self['py']).nn_interpolator(zz)\
                 [LE[0]:RE[0]:side*1j, \
                  LE[1]:RE[1]:side*1j]
        if log_spacing:
            zi = 10**(zi)
        return [xi,yi,zi]

    _okay_to_serialize = True

    def _store_fields(self, fields, node_name = None, force = False):
        fields = ensure_list(fields)
        if node_name is None: node_name = self._gen_node_name()
        for field in fields:
            #mylog.debug("Storing %s in node %s",
                #self._convert_field_name(field), node_name)
            self.hierarchy.save_data(self[field], node_name,
                self._convert_field_name(field), force = force,
                passthrough = True)

    def _obtain_fields(self, fields, node_name = None):
        if not self._okay_to_serialize: return
        fields = ensure_list(fields)
        if node_name is None: node_name = self._gen_node_name()
        for field in fields:
            #mylog.debug("Trying to obtain %s from node %s",
                #self._convert_field_name(field), node_name)
            fdata=self.hierarchy.get_data(node_name, 
                self._convert_field_name(field))
            if fdata is not None:
                #mylog.debug("Got %s from node %s", field, node_name)
                self[field] = fdata[:]
        return True

    def _deserialize(self, node_name = None):
        if not self._okay_to_serialize: return
        self._obtain_fields(self._key_fields, node_name)
        self._obtain_fields(self.fields, node_name)

    def _serialize(self, node_name = None, force = False):
        if not self._okay_to_serialize: return
        self._store_fields(self._key_fields, node_name, force)
        self._store_fields(self.fields, node_name, force)

class AMRSliceBase(AMR2DData):
    """
    AMRSlice is an orthogonal slice through the data, taking all the points
    at the finest resolution available and then indexing them.  It is more
    appropriately thought of as a slice 'operator' than an object,
    however, as its field and coordinate can both change.
    """

    _top_node = "/Slices"
    _type_name = "slice"
    _con_args = ('axis', 'coord')
    #@time_execution
    def __init__(self, axis, coord, fields = None, center=None, pf=None,
                 node_name = False, **kwargs):
        """
        Slice along *axis*:ref:`axis-specification`, at the coordinate *coord*.
        Optionally supply fields.
        """
        AMR2DData.__init__(self, axis, fields, pf, **kwargs)
        self.center = center
        if center is not None: self.set_field_parameter('center',center)
        self.coord = coord
        if node_name is False:
            self._refresh_data()
        else:
            if node_name is True: self._deserialize()
            else: self._deserialize(node_name)

    def reslice(self, coord):
        """
        Change the entire dataset, clearing out the current data and slicing at
        a new location.  Not terribly useful except for in-place plot changes.
        """
        mylog.debug("Setting coordinate to %0.5e" % coord)
        self.coord = coord
        self._refresh_data()

    def shift(self, val):
        """
        Moves the slice coordinate up by either a floating point value, or an
        integer number of indices of the finest grid.
        """
        if isinstance(val, types.FloatType):
            # We add the dx
            self.coord += val
        elif isinstance(val, types.IntType):
            # Here we assume that the grid is the max level
            level = self.hierarchy.max_level
            self.coord
            dx = self.hierarchy.gridDxs[self.hierarchy.levelIndices[level][0]]
            self.coord += dx * val
        else:
            raise ValueError(val)
        self._refresh_data()

    def _generate_coords(self):
        points = []
        for grid in self._get_grids():
            points.append(self._generate_grid_coords(grid))
        t = self._mpi_catarray(na.concatenate(points))
        self['px'] = t[:,0]
        self['py'] = t[:,1]
        self['pz'] = t[:,2]
        self['pdx'] = t[:,3]
        self['pdy'] = t[:,4]
        self['pdz'] = t[:,3] # Does not matter!

        # Now we set the *actual* coordinates
        self[axis_names[x_dict[self.axis]]] = t[:,0]
        self[axis_names[y_dict[self.axis]]] = t[:,1]
        self[axis_names[self.axis]] = t[:,2]

        self.ActiveDimensions = (t.shape[0], 1, 1)

    def _get_list_of_grids(self):
        goodI = ((self.pf.h.gridRightEdge[:,self.axis] > self.coord)
              &  (self.pf.h.gridLeftEdge[:,self.axis] <= self.coord ))
        self._grids = self.pf.h.grids[goodI] # Using sources not hierarchy

    def __cut_mask_child_mask(self, grid):
        mask = grid.child_mask.copy()
        return mask

    def _generate_grid_coords(self, grid):
        xaxis = x_dict[self.axis]
        yaxis = y_dict[self.axis]
        ds, dx, dy = grid.dds[self.axis], grid.dds[xaxis], grid.dds[yaxis]
        wantedIndex = int(((self.coord-grid.LeftEdge[self.axis])/ds))
        sl = [slice(None), slice(None), slice(None)]
        sl[self.axis] = slice(wantedIndex, wantedIndex + 1)
        #sl.reverse()
        sl = tuple(sl)
        nx = grid.child_mask.shape[xaxis]
        ny = grid.child_mask.shape[yaxis]
        mask = self.__cut_mask_child_mask(grid)[sl]
        cm = na.where(mask.ravel()== 1)
        cmI = na.indices((nx,ny))
        xind = cmI[0,:].ravel()
        xpoints = na.ones(cm[0].shape, 'float64')
        xpoints *= xind[cm]*dx+(grid.LeftEdge[xaxis] + 0.5*dx)
        yind = cmI[1,:].ravel()
        ypoints = na.ones(cm[0].shape, 'float64')
        ypoints *= yind[cm]*dy+(grid.LeftEdge[yaxis] + 0.5*dy)
        zpoints = na.ones(xpoints.shape, 'float64') * self.coord
        dx = na.ones(xpoints.shape, 'float64') * dx/2.0
        dy = na.ones(xpoints.shape, 'float64') * dy/2.0
        t = na.array([xpoints, ypoints, zpoints, dx, dy]).swapaxes(0,1)
        return t

    @restore_grid_state
    def _get_data_from_grid(self, grid, field):
        # So what's our index of slicing?  This is what we need to figure out
        # first, so we can deal with our data in the fastest way.
        dx = grid.dds[self.axis]
        wantedIndex = int(((self.coord-grid.LeftEdge[self.axis])/dx))
        sl = [slice(None), slice(None), slice(None)]
        sl[self.axis] = slice(wantedIndex, wantedIndex + 1)
        sl = tuple(sl)
        if self.pf.field_info.has_key(field) and self.pf.field_info[field].particle_type:
            return grid[field]
        elif field in self.pf.field_info and self.pf.field_info[field].not_in_all:
            dv = grid[field][sl]
        elif not grid.has_key(field):
            conv_factor = 1.0
            if self.pf.field_info.has_key(field):
                conv_factor = self.pf.field_info[field]._convert_function(self)
            dv = self._read_data_slice(grid, field, self.axis, wantedIndex) * conv_factor
        else:
            dv = grid[field]
            if dv.size == 1: dv = na.ones(grid.ActiveDimensions)*dv
            dv = dv[sl]
        mask = self.__cut_mask_child_mask(grid)[sl]
        dataVals = dv.ravel()[mask.ravel() == 1]
        return dataVals

    def _gen_node_name(self):
        return "%s/%s_%s" % \
            (self._top_node, self.axis, self.coord)

class AMRCuttingPlaneBase(AMR2DData):
    """
    AMRCuttingPlane is an oblique plane through the data,
    defined by a normal vector and a coordinate.  It attempts to guess
    an 'up' vector, which cannot be overridden, and then it pixelizes
    the appropriate data onto the plane without interpolation.
    """
    _plane = None
    _top_node = "/CuttingPlanes"
    _key_fields = AMR2DData._key_fields + ['pz','pdz']
    _type_name = "cutting"
    _con_args = ('normal', 'center')
    def __init__(self, normal, center, fields = None, node_name = None,
                 **kwargs):
        """
        The Cutting Plane slices at an oblique angle, where we use
        the *normal* vector and the *center* to define the viewing plane.
        The 'up' direction is guessed at automatically.
        """
        AMR2DData.__init__(self, 4, fields, **kwargs)
        self.center = center
        self.set_field_parameter('center',center)
        # Let's set up our plane equation
        # ax + by + cz + d = 0
        self._norm_vec = normal/na.sqrt(na.dot(normal,normal))
        self._d = -1.0 * na.dot(self._norm_vec, self.center)
        # First we try all three, see which has the best result:
        vecs = na.identity(3)
        _t = na.cross(self._norm_vec, vecs).sum(axis=1)
        ax = _t.argmax()
        self._x_vec = na.cross(vecs[ax,:], self._norm_vec).ravel()
        self._x_vec /= na.sqrt(na.dot(self._x_vec, self._x_vec))
        self._y_vec = na.cross(self._norm_vec, self._x_vec).ravel()
        self._y_vec /= na.sqrt(na.dot(self._y_vec, self._y_vec))
        self._rot_mat = na.array([self._x_vec,self._y_vec,self._norm_vec])
        self._inv_mat = na.linalg.pinv(self._rot_mat)
        self.set_field_parameter('cp_x_vec',self._x_vec)
        self.set_field_parameter('cp_y_vec',self._y_vec)
        self.set_field_parameter('cp_z_vec',self._norm_vec)
        if node_name is False:
            self._refresh_data()
        else:
            if node_name is True: self._deserialize()
            else: self._deserialize(node_name)

    @property
    def normal(self):
        return self._norm_vec

    def _get_list_of_grids(self):
        # Recall that the projection of the distance vector from a point
        # onto the normal vector of a plane is:
        # D = (a x_0 + b y_0 + c z_0 + d)/sqrt(a^2+b^2+c^2)
        # @todo: Convert to using corners
        LE = self.pf.h.gridLeftEdge
        RE = self.pf.h.gridRightEdge
        vertices = na.array([[LE[:,0],LE[:,1],LE[:,2]],
                             [RE[:,0],RE[:,1],RE[:,2]],
                             [LE[:,0],LE[:,1],RE[:,2]],
                             [RE[:,0],RE[:,1],LE[:,2]],
                             [LE[:,0],RE[:,1],RE[:,2]],
                             [RE[:,0],LE[:,1],LE[:,2]],
                             [LE[:,0],RE[:,1],LE[:,2]],
                             [RE[:,0],LE[:,1],RE[:,2]]])
        # This gives us shape: 8, 3, n_grid
        D = na.sum(self._norm_vec.reshape((1,3,1)) * vertices, axis=1) + self._d
        self.D = D
        self._grids = self.hierarchy.grids[
            na.where(na.logical_not(na.all(D<0,axis=0) | na.all(D>0,axis=0) )) ]

    @cache_mask
    def _get_cut_mask(self, grid):
        # This is slow.  Suggestions for improvement would be great...
        ss = grid.ActiveDimensions
        D = na.ones(ss) * self._d
        x = grid.LeftEdge[0] + grid.dds[0] * \
                (na.arange(grid.ActiveDimensions[0], dtype='float64')+0.5)
        y = grid.LeftEdge[1] + grid.dds[1] * \
                (na.arange(grid.ActiveDimensions[1], dtype='float64')+0.5)
        z = grid.LeftEdge[2] + grid.dds[2] * \
                (na.arange(grid.ActiveDimensions[2], dtype='float64')+0.5)
        D += (x * self._norm_vec[0]).reshape(ss[0],1,1)
        D += (y * self._norm_vec[1]).reshape(1,ss[1],1)
        D += (z * self._norm_vec[2]).reshape(1,1,ss[2])
        diag_dist = na.sqrt(na.sum(grid.dds**2.0))
        cm = (na.abs(D) <= 0.5*diag_dist) # Boolean
        return cm

    def _generate_coords(self):
        points = []
        for grid in self._get_grids():
            points.append(self._generate_grid_coords(grid))
        t = self._mpi_catarray(na.concatenate(points))
        pos = (t[:,0:3] - self.center)
        self['px'] = na.dot(pos, self._x_vec)
        self['py'] = na.dot(pos, self._y_vec)
        self['pz'] = na.dot(pos, self._norm_vec)
        self['pdx'] = t[:,3] * 0.5
        self['pdy'] = t[:,3] * 0.5
        self['pdz'] = t[:,3] * 0.5

    def _generate_grid_coords(self, grid):
        pointI = self._get_point_indices(grid)
        coords = [grid[ax][pointI].ravel() for ax in 'xyz']
        coords.append(na.ones(coords[0].shape, 'float64') * just_one(grid['dx']))
        return na.array(coords).swapaxes(0,1)

    def _get_data_from_grid(self, grid, field):
        if not self.pf.field_info[field].particle_type:
            pointI = self._get_point_indices(grid)
            if grid[field].size == 1: # dx, dy, dz, cellvolume
                t = grid[field] * na.ones(grid.ActiveDimensions)
                return t[pointI].ravel()
            return grid[field][pointI].ravel()
        else:
            return grid[field]

    def interpolate_discretize(self, *args, **kwargs):
        pass

    @cache_point_indices
    def _get_point_indices(self, grid, use_child_mask=True):
        k = na.zeros(grid.ActiveDimensions, dtype='bool')
        k = (k | self._get_cut_mask(grid))
        if use_child_mask: k = (k & grid.child_mask)
        return na.where(k)

    def _gen_node_name(self):
        cen_name = ("%s" % self.center).replace(" ","_")[1:-1]
        L_name = ("%s" % self._norm_vec).replace(" ","_")[1:-1]
        return "%s/c%s_L%s" % \
            (self._top_node, cen_name, L_name)

class AMRProjBase(AMR2DData):
    _top_node = "/Projections"
    _key_fields = AMR2DData._key_fields + ['weight_field']
    _type_name = "proj"
    _con_args = ('axis', 'field', 'weight_field')
    def __init__(self, axis, field, weight_field = None,
                 max_level = None, center = None, pf = None,
                 source=None, node_name = None, field_cuts = None, **kwargs):
        """
        AMRProj is a projection of a *field* along an *axis*.  The field
        can have an associated *weight_field*, in which case the values are
        multiplied by a weight before being summed, and then divided by the sum
        of that weight.
        """
        AMR2DData.__init__(self, axis, field, pf, node_name = None, **kwargs)
        self._field_cuts = field_cuts
        self.center = center
        if center is not None: self.set_field_parameter('center',center)
        self._node_name = node_name
        self._initialize_source(source)
        self._grids = self.source._grids
        if max_level == None:
            max_level = self.hierarchy.maxLevel
        if self.source is not None:
            max_level = min(max_level, self.source.gridLevels.max())
        self._max_level = max_level
        self._weight = weight_field
        self.func = na.sum # for the future
        self.__retval_coords = {}
        self.__retval_fields = {}
        self.__retval_coarse = {}
        self.__overlap_masks = {}
        self._temp = {}
        self._deserialize(node_name)
        self._refresh_data()
        if self._okay_to_serialize: self._serialize(node_name=self._node_name)

    def _convert_field_name(self, field):
        if field == "weight_field": return "weight_field_%s" % self._weight
        if field in self._key_fields: return field
        return "%s_%s" % (field, self._weight)

    def _initialize_source(self, source = None):
        if source is None:
            check, source = self._partition_hierarchy_2d(self.axis)
            self._check_region = check
            #self._okay_to_serialize = (not check)
        else:
            self._distributed = False
            self._okay_to_serialize = False
            self._check_region = True
        self.source = source
        if self._field_cuts is not None:
            # Override if field cuts are around; we don't want to serialize!
            self._check_region = True
            self._okay_to_serialize = False
        if self._node_name is not None:
            self._node_name = "%s/%s" % (self._top_node,self._node_name)
            self._okay_to_serialize = True

    #@time_execution
    def __cache_data(self):
        rdf = self.hierarchy.grid.readDataFast
        self.hierarchy.grid.readDataFast = readDataPackedHandle
        for fn,g_list in self.hierarchy.cpu_map.items():
            to_read = na.intersect1d(g_list, self.source._grids)
            if len(to_read) == 0: continue
            fh = tables.openFile(to_read[0].filename,'r')
            for g in to_read:
                g.handle = fh
                for field in ensure_list(self.fields):
                    g[field]
                del g.handle
            fh.close()
        self.hierarchy.grid.readDataFast = readDataPackedHandle

    #@time_execution
    def __calculate_overlap(self):
        s = self.source
        mylog.info("Generating overlap masks")
        i = 0
        pbar = get_pbar("Reading and masking grids ", len(s._grids))
        for level in range(self._max_level+1):
            mylog.debug("Examining level %s", level)
            grids = s.levelIndices[level]
            RE = s.gridRightEdge[grids]
            LE = s.gridLeftEdge[grids]
            for grid in s._grids[grids]:
                pbar.update(i)
                self.__overlap_masks[grid.id] = \
                    grid._generate_overlap_masks(self.axis, LE, RE)
                i += 1
        pbar.finish()
        mylog.info("Finished calculating overlap.")

    def __get_dls(self, grid, fields):
        # Place holder for a time when maybe we will not be doing just
        # a single dx for every field.
        dls = []
        convs = []
        for field in fields + [self._weight]:
            if field is None: continue
            dls.append(just_one(grid['d%s' % axis_names[self.axis]]))
            convs.append(self.pf.units[self.pf.field_info[field].projection_conversion])
        return na.array(dls), na.array(convs)

    def __project_level(self, level, fields):
        grids_to_project = self.source.select_grids(level)
        dls, convs = self.__get_dls(grids_to_project[0], fields)
        zero_out = (level != self._max_level)
        pbar = get_pbar('Projecting  level % 2i / % 2i ' \
                          % (level, self._max_level), len(grids_to_project))
        for pi, grid in enumerate(grids_to_project):
            g_coords, g_fields = self._project_grid(grid, fields, zero_out)
            self.__retval_coords[grid.id] = g_coords
            self.__retval_fields[grid.id] = g_fields
            for fi in range(len(fields)): g_fields[fi] *= dls[fi]
            if self._weight is not None: g_coords[3] *= dls[-1]
            pbar.update(pi)
            grid.clear_data()
        pbar.finish()
        self.__combine_grids_on_level(level) # In-place
        if level > 0 and level <= self._max_level:
            self.__refine_to_level(level) # In-place
        coord_data = []
        field_data = []
        for grid in grids_to_project:
            coarse = self.__retval_coords[grid.id][2]==0 # Where childmask = 0
            fine = ~coarse
            coord_data.append([pi[fine] for pi in self.__retval_coords[grid.id]])
            field_data.append([pi[fine] for pi in self.__retval_fields[grid.id]])
            self.__retval_coords[grid.id] = [pi[coarse] for pi in self.__retval_coords[grid.id]]
            self.__retval_fields[grid.id] = [pi[coarse] for pi in self.__retval_fields[grid.id]]
        coord_data = na.concatenate(coord_data, axis=1)
        field_data = na.concatenate(field_data, axis=1)
        if self._weight is not None:
            field_data = field_data / coord_data[3,:].reshape((1,coord_data.shape[1]))
        else:
            field_data *= convs[...,na.newaxis]
        mylog.info("Level %s done: %s final", \
                   level, coord_data.shape[1])
        dx = grids_to_project[0].dds[self.axis] # this is our dl
        return coord_data, dx, field_data

    def __combine_grids_on_level(self, level):
        grids = self.source.select_grids(level)
        grids_i = self.source.levelIndices[level]
        pbar = get_pbar('Combining   level % 2i / % 2i ' \
                          % (level, self._max_level), len(grids))
        # We have an N^2 check, so we try to be as quick as possible
        # and to skip as many as possible
        for pi, grid1 in enumerate(grids):
            pbar.update(pi)
            if self.__retval_coords[grid1.id][0].shape[0] == 0: continue
            for grid2 in self.source._grids[grids_i][self.__overlap_masks[grid1.id]]:
                if self.__retval_coords[grid2.id][0].shape[0] == 0 \
                  or grid1.id == grid2.id:
                    continue
                args = [] # First is source, then destination
                args += self.__retval_coords[grid2.id] + [self.__retval_fields[grid2.id]]
                args += self.__retval_coords[grid1.id] + [self.__retval_fields[grid1.id]]
                args.append(1) # Refinement factor
                args.append(na.ones(args[0].shape, dtype='int64'))
                kk = PointCombine.CombineGrids(*args)
                goodI = args[-1].astype('bool')
                self.__retval_coords[grid2.id] = \
                    [coords[goodI] for coords in self.__retval_coords[grid2.id]]
                self.__retval_fields[grid2.id] = \
                    [fields[goodI] for fields in self.__retval_fields[grid2.id]]
        pbar.finish()

    def __refine_to_level(self, level):
        grids = self.source.select_grids(level)
        grids_up = self.source.levelIndices[level-1]
        pbar = get_pbar('Refining to level % 2i / % 2i ' \
                          % (level, self._max_level), len(grids))
        for pi, grid1 in enumerate(grids):
            pbar.update(pi)
            for parent in ensure_list(grid1.Parent):
                if parent.id not in self.__overlap_masks: continue
                for grid2 in self.source._grids[grids_up][self.__overlap_masks[parent.id]]:
                    if self.__retval_coords[grid2.id][0].shape[0] == 0: continue
                    args = []
                    args += self.__retval_coords[grid2.id] + [self.__retval_fields[grid2.id]]
                    args += self.__retval_coords[grid1.id] + [self.__retval_fields[grid1.id]]
                    # Refinement factor, which is same in all directions
                    args.append(int(grid2['dx'] / grid1['dx'])) 
                    args.append(na.ones(args[0].shape, dtype='int64'))
                    kk = PointCombine.CombineGrids(*args)
                    goodI = args[-1].astype('bool')
                    self.__retval_coords[grid2.id] = \
                        [coords[goodI] for coords in self.__retval_coords[grid2.id]]
                    self.__retval_fields[grid2.id] = \
                        [fields[goodI] for fields in self.__retval_fields[grid2.id]]
        for grid1 in self.source.select_grids(level-1):
            if not self._check_region and self.__retval_coords[grid1.id][0].size != 0:
                mylog.error("Something messed up, and %s still has %s points of data",
                            grid1, self.__retval_coords[grid1.id][0].size)
                mylog.error("You might try setting the ReconstructHierarchy option in [lagos]")
                raise ValueError(grid1, self.__retval_coords[grid1.id])
        pbar.finish()

    #@time_execution
    def get_data(self, fields = None):
        if fields is None: fields = ensure_list(self.fields)[:]
        else: fields = ensure_list(fields)
        self._obtain_fields(fields, self._node_name)
        fields = [f for f in fields if f not in self.data]
        if len(fields) == 0: return
        self.__calculate_overlap()
        coord_data = []
        field_data = []
        dxs = []
        # We do this here, but I am not convinced it should be done here
        # It is probably faster, as it consolidates IO, but if we did it in
        # _project_level, then it would be more memory conservative
        self._preload(self.source._grids, self._get_dependencies(fields),
                      self.hierarchy.queue)
        for level in range(0, self._max_level+1):
            my_coords, my_dx, my_fields = self.__project_level(level, fields)
            coord_data.append(my_coords)
            field_data.append(my_fields)
            dxs.append(my_dx * na.ones(my_coords.shape[1], dtype='float64'))
            if self._check_region and False:
                check=self.__cleanup_level(level - 1)
                if len(check) > 0: all_data.append(check)
            # Now, we should clean up after ourselves...
            for grid in self.source.select_grids(level - 1):
                grid.clear_data()
                del self.__retval_coords[grid.id]
                del self.__retval_fields[grid.id]
                del self.__overlap_masks[grid.id]
        coord_data = na.concatenate(coord_data, axis=1)
        field_data = na.concatenate(field_data, axis=1)
        dxs = na.concatenate(dxs, axis=1)
        # We now convert to half-widths and center-points
        data = {}
        data['pdx'] = dxs
        data['px'] = (coord_data[0,:]+0.5) * data['pdx']
        data['py'] = (coord_data[1,:]+0.5) * data['pdx']
        data['weight_field'] = coord_data[3,:].copy()
        del coord_data
        data['pdx'] *= 0.5
        data['pdy'] = data['pdx'] # generalization is out the window!
        data['fields'] = field_data
        # Now we run the finalizer, which is ignored if we don't need it
        data = self._mpi_catdict(data)
        field_data = na.vsplit(data.pop('fields'), len(fields))
        for fi, field in enumerate(fields):
            self[field] = field_data[fi].ravel()
            self._store_fields(field, self._node_name)
        for i in data.keys(): self[i] = data.pop(i)

    def add_fields(self, fields, weight = "CellMassMsun"):
        pass

    def _project_grid(self, grid, fields, zero_out):
        if self._weight is None:
            weight_data = na.ones(grid.ActiveDimensions, dtype='float64')
        else:
            weight_data = self._get_data_from_grid(grid, self._weight).astype('float64')
        if zero_out: weight_data[grid.child_indices] = 0
        # if we zero it out here, then we only have to zero out the weight!
        masked_data = [self._get_data_from_grid(grid, field) * weight_data
                       for field in fields]
        full_proj = [self.func(field,axis=self.axis) for field in masked_data]
        weight_proj = self.func(weight_data,axis=self.axis)
        if (self._check_region and not self.source._is_fully_enclosed(grid)) or self._field_cuts is not None:
            used_data = self._get_points_in_region(grid).astype('bool')
            used_points = na.where(na.logical_or.reduce(used_data, self.axis))
        else:
            used_data = na.array([1.0], dtype='bool')
            used_points = slice(None)
        if zero_out:
            subgrid_mask = na.logical_and.reduce(
                                na.logical_or(grid.child_mask,
                                             ~used_data),
                                self.axis).astype('int64')
        else:
            subgrid_mask = na.ones(full_proj[0].shape, dtype='int64')
        xind, yind = [arr[used_points].ravel() for arr in na.indices(full_proj[0].shape)]
        start_index = grid.get_global_startindex()
        xpoints = (xind + (start_index[x_dict[self.axis]])).astype('int64')
        ypoints = (yind + (start_index[y_dict[self.axis]])).astype('int64')
        return ([xpoints, ypoints,
                subgrid_mask[used_points].ravel(),
                weight_proj[used_points].ravel()],
                [data[used_points].ravel() for data in full_proj])

    def _get_points_in_region(self, grid):
        pointI = self.source._get_point_indices(grid, use_child_mask=False)
        point_mask = na.zeros(grid.ActiveDimensions)
        point_mask[pointI] = 1.0
        if self._field_cuts is not None:
            for cut in self._field_cuts:
                point_mask *= eval(cut)
        return point_mask

    @restore_grid_state
    def _get_data_from_grid(self, grid, field):
        if self._check_region:
            bad_points = self._get_points_in_region(grid)
        else:
            bad_points = 1.0
        d = grid[field] * bad_points
        if grid.id == 1: self._temp[grid.id] = d
        return d

    def _gen_node_name(self):
        return  "%s/%s" % \
            (self._top_node, self.axis)

class AMR3DData(AMRData, GridPropertiesMixin):
    _key_fields = ['x','y','z','dx','dy','dz']
    """
    Class describing a cluster of data points, not necessarily sharing any
    particular attribute.
    """
    _spatial = False
    _num_ghost_zones = 0
    def __init__(self, center, fields, pf = None, **kwargs):
        """
        Returns an instance of AMR3DData, or prepares one.  Usually only
        used as a base class.  Note that *center* is supplied, but only used
        for fields and quantities that require it.
        """
        AMRData.__init__(self, pf, fields, **kwargs)
        self.center = center
        self.set_field_parameter("center",center)
        self.coords = None
        self._grids = None

    def _generate_coords(self):
        mylog.info("Generating coords for %s grids", len(self._grids))
        points = []
        for i,grid in enumerate(self._grids):
            #grid._generate_coords()
            if ( (i%100) == 0):
                mylog.info("Working on % 7i / % 7i", i, len(self._grids))
            grid.set_field_parameter("center", self.center)
            points.append((na.ones(
                grid.ActiveDimensions,dtype='float64')*grid['dx'])\
                    [self._get_point_indices(grid)])
            t = na.concatenate([t,points])
            del points
        self['dx'] = t
        #self['dy'] = t
        #self['dz'] = t
        mylog.info("Done with coordinates")

    @restore_grid_state
    def _generate_grid_coords(self, grid, field=None):
        pointI = self._get_point_indices(grid)
        dx = na.ones(pointI[0].shape[0], 'float64') * grid.dds[0]
        tr = na.array([grid['x'][pointI].ravel(), \
                grid['y'][pointI].ravel(), \
                grid['z'][pointI].ravel(), \
                grid["RadiusCode"][pointI].ravel(),
                dx, grid["GridIndices"][pointI].ravel()], 'float64').swapaxes(0,1)
        return tr

    def get_data(self, fields=None, in_grids=False):
        if self._grids == None:
            self._get_list_of_grids()
        points = []
        #if not self.has_key('dx'):
            #self._generate_coords()
        if not fields:
            fields_to_get = self.fields
        else:
            fields_to_get = ensure_list(fields)
        mylog.debug("Going to obtain %s", fields_to_get)
        for field in fields_to_get:
            if self.data.has_key(field):
                continue
            mylog.info("Getting field %s from %s", field, len(self._grids))
            if field not in self.hierarchy.field_list and not in_grids:
                if self._generate_field(field):
                    continue # True means we already assigned it
            self[field] = na.concatenate(
                [self._get_data_from_grid(grid, field)
                 for grid in self._grids])

    @restore_grid_state
    def _get_data_from_grid(self, grid, field):
        if field in self.pf.field_info and self.pf.field_info[field].particle_type:
            # int64 -> float64 with the first real set of data
            if grid.NumberOfParticles == 0: return na.array([], dtype='int64')
            pointI = self._get_particle_indices(grid)
            if self.pf.field_info[field].vector_field:
                f = grid[field]
                return na.array([f[i,:][pointI] for i in range(3)])
            return grid[field][pointI].ravel()
        if field in self.pf.field_info and self.pf.field_info[field].vector_field:
            pointI = self._get_point_indices(grid)
            f = grid[field]
            return na.array([f[i,:][pointI] for i in range(3)])
        else:
            pointI = self._get_point_indices(grid)
            if grid[field].size == 1: # dx, dy, dz, cellvolume
                t = grid[field] * na.ones(grid.ActiveDimensions, dtype='float64')
                return t[pointI].ravel()
            return grid[field][pointI].ravel()

    def _flush_data_to_grids(self, field, default_val, dtype='float32'):
        """
        A dangerous, thusly underscored, thing to do to a data object,
        we can flush back any changes in a given field that have been made
        with a default value for the rest of the grid.
        """
        i = 0
        for grid in self._grids:
            pointI = self._get_point_indices(grid)
            np = pointI[0].ravel().size
            new_field = na.ones(grid.ActiveDimensions, dtype=dtype) * default_val
            new_field[pointI] = self[field][i:i+np]
            if grid.data.has_key(field): del grid.data[field]
            grid[field] = new_field
            i += np

    def _generate_field(self, field):
        if self.pf.field_info.has_key(field):
            # First we check the validator
            try:
                self.pf.field_info[field].check_available(self)
            except NeedsGridType, ngt_exception:
                # We leave this to be implementation-specific
                self._generate_field_in_grids(field, ngt_exception.ghost_zones)
                return False
            else:
                self[field] = self.pf.field_info[field](self)
                return True
        else: # Can't find the field, try as it might
            raise exceptions.KeyError(field)

    def _generate_field_in_grids(self, field, num_ghost_zones=0):
        for grid in self._grids:
            self.__touch_grid_field(grid, field)

    @restore_grid_state
    def __touch_grid_field(self, grid, field):
        grid[field]

    def _is_fully_enclosed(self, grid):
        return na.all(self._get_cut_mask)

    def _get_point_indices(self, grid, use_child_mask=True):
        k = na.zeros(grid.ActiveDimensions, dtype='bool')
        k = (k | self._get_cut_mask(grid))
        if use_child_mask: k = (k & grid.child_mask)
        return na.where(k)

    def _get_cut_particle_mask(self, grid):
        fake_grid = FakeGridForParticles(grid)
        return self._get_cut_mask(fake_grid)

    def _get_particle_indices(self, grid):
        k = na.zeros(grid.NumberOfParticles, dtype='bool')
        k = (k | self._get_cut_particle_mask(grid))
        return na.where(k)

    def cut_region(self, field_cuts):
        """
        Return an InLineExtractedRegion, where the grid cells are cut on the
        fly with a set of field_cuts.
        """
        return InLineExtractedRegionBase(self, field_cuts)

    def extract_region(self, indices):
        """
        Return an ExtractedRegion where the points contained in it are defined
        as the points in `this` data object with the given *indices*.
        """
        return ExtractedRegionBase(self, indices)

    def __get_quantities(self):
        if self.__quantities is None:
            self.__quantities = DerivedQuantityCollection(self)
        return self.__quantities
    __quantities = None
    quantities = property(__get_quantities)

    def extract_connected_sets(self, field, num_levels, min_val, max_val,
                                log_space=True, cumulative=True, cache=False):
        """
        This function will create a set of contour objects, defined
        by having connected cell structures, which can then be
        studied and used to 'paint' their source grids, thus enabling
        them to be plotted.
        """
        if log_space:
            cons = na.logspace(na.log10(min_val),na.log10(max_val),
                               num_levels+1)
        else:
            cons = na.linspace(min_val, max_val, num_levels+1)
        contours = {}
        if cache: cached_fields = defaultdict(lambda: dict())
        else: cached_fields = None
        for level in range(num_levels):
            contours[level] = {}
            if cumulative:
                mv = max_val
            else:
                mv = cons[level+1]
            cids = identify_contours(self, field, cons[level], mv,
                                     cached_fields)
            for cid, cid_ind in cids.items():
                contours[level][cid] = self.extract_region(cid_ind)
        return cons, contours

    def paint_grids(self, field, value, default_value=None):
        """
        This function paints every cell in our dataset with a given *value*.
        If default_value is given, the other values for the given in every grid
        are discarded and replaced with *default_value*.  Otherwise, the field is
        mandated to 'know how to exist' in the grid.

        Note that this only paints the cells *in the dataset*, so cells in grids
        with child cells are left untouched.
        """
        for grid in self._grids:
            if default_value != None:
                grid[field] = na.ones(grid.ActiveDimensions)*value
            grid[field][self._get_point_indices()] = value


class ExtractedRegionBase(AMR3DData):
    """
    ExtractedRegions are arbitrarily defined containers of data, useful
    for things like selection along a baryon field.
    """
    _type_name = "extracted_region"
    _con_args = ('_base_region', '_indices')
    def __init__(self, base_region, indices, force_refresh=True, **kwargs):
        cen = base_region.get_field_parameter("center")
        AMR3DData.__init__(self, center=cen,
                            fields=None, pf=base_region.pf, **kwargs)
        self._base_region = base_region # We don't weakly reference because
                                        # It is not cyclic
        if isinstance(indices, types.DictType):
            self._indices = indices
            self._grids = self._base_region.pf.h.grids[self._indices.keys()]
        else:
            self._grids = None
            self._base_indices = indices
        if force_refresh: self._refresh_data()

    def _get_cut_particle_mask(self, grid):
        # Override to provide a warning
        mylog.warning("Returning all particles from an Extracted Region.  This could be incorrect!")
        return True

    def _get_list_of_grids(self):
        # Okay, so what we're going to want to do is get the pointI from
        # region._get_point_indices(grid) for grid in base_region._grids,
        # and then construct an array of those, which we will select along indices.
        if self._grids != None: return
        grid_vals, xi, yi, zi = [], [], [], []
        for grid in self._base_region._grids:
            xit,yit,zit = self._base_region._get_point_indices(grid)
            grid_vals.append(na.ones(xit.shape, dtype='int') * (grid.id-grid._id_offset))
            xi.append(xit)
            yi.append(yit)
            zi.append(zit)
        grid_vals = na.concatenate(grid_vals)[self._base_indices]
        grid_order = na.argsort(grid_vals)
        # Note: grid_vals is still unordered
        grid_ids = na.unique(grid_vals)
        xi = na.concatenate(xi)[self._base_indices][grid_order]
        yi = na.concatenate(yi)[self._base_indices][grid_order]
        zi = na.concatenate(zi)[self._base_indices][grid_order]
        bc = na.bincount(grid_vals)
        splits = []
        for i,v in enumerate(bc):
            if v > 0: splits.append(v)
        splits = na.add.accumulate(splits)
        xis, yis, zis = [na.array_split(aa, splits) for aa in [xi,yi,zi]]
        self._indices = {}
        h = self._base_region.pf.h
        for grid_id, x, y, z in zip(grid_ids, xis, yis, zis):
            # grid_id needs no offset
            ll = h.grids[grid_id].ActiveDimensions.prod() \
               - (na.logical_not(h.grids[grid_id].child_mask)).sum()
            # This means we're completely enclosed, except for child masks
            if x.size == ll:
                self._indices[grid_id] = None
            else:
                # This will slow things down a bit, but conserve memory
                self._indices[grid_id] = \
                    na.zeros(h.grids[grid_id].ActiveDimensions, dtype='bool')
                self._indices[grid_id][(x,y,z)] = True
        self._grids = h.grids[self._indices.keys()]

    def _is_fully_enclosed(self, grid):
        if self._indices[grid.id-grid._id_offset] is None or \
            (self._indices[grid.id-grid._id_offset][0].size ==
             grid.ActiveDimensions.prod()):
            return True
        return False

    __empty_array = na.array([], dtype='bool')
    def _get_point_indices(self, grid, use_child_mask=True):
        # Yeah, if it's not true, we don't care.
        tr = self._indices.get(grid.id-grid._id_offset, self.__empty_array)
        if tr is None: tr = na.where(grid.child_mask)
        else: tr = na.where(tr)
        return tr

    def __repr__(self):
        # We'll do this the slow way to be clear what's going on
        s = "%s (%s): " % (self.__class__.__name__, self.pf)
        s += ", ".join(["%s=%s" % (i, getattr(self,i))
                       for i in self._con_args if i != "_indices"])
        return s

    def join(self, other):
        ng = {}
        gs = set(self._indices.keys() + other._indices.keys())
        for g in gs:
            grid = self.pf.h.grids[g]
            if g in other._indices and g in self._indices:
                # We now join the indices
                ind = na.zeros(grid.ActiveDimensions, dtype='bool')
                ind[self._indices[g]] = True
                ind[other._indices[g]] = True
                if ind.prod() == grid.ActiveDimensions.prod(): ind = None
            elif g in self._indices:
                ind = self._indices[g]
            elif g in other._indices:
                ind = other._indices[g]
            # Okay we have indices
            if ind is not None: ind = ind.copy()
            ng[g] = ind
        gl = self.pf.h.grids[list(gs)]
        gc = self.pf.h.grid_collection(
            self._base_region.get_field_parameter("center"), gl)
        return self.pf.h.extracted_region(gc, ng)

class InLineExtractedRegionBase(AMR3DData):
    """
    In-line extracted regions accept a base region and a set of field_cuts to
    determine which points in a grid should be included.
    """
    def __init__(self, base_region, field_cuts, **kwargs):
        cen = base_region.get_field_parameter("center")
        AMR3DData.__init__(self, center=cen,
                            fields=None, pf=base_region.pf, **kwargs)
        self._base_region = base_region # We don't weakly reference because
                                        # It is not cyclic
        self._field_cuts = ensure_list(field_cuts)[:]
        self._refresh_data()

    def _get_list_of_grids(self):
        self._grids = self._base_region._grids

    def _is_fully_enclosed(self, grid):
        return False

    @cache_mask
    def _get_cut_mask(self, grid):
        point_mask = na.ones(grid.ActiveDimensions, dtype='bool')
        point_mask *= self._base_region._get_cut_mask(grid)
        for cut in self._field_cuts:
            point_mask *= eval(cut)
        return point_mask

class AMRCylinderBase(AMR3DData):
    """
    We can define a cylinder (or disk) to act as a data object.
    """
    _type_name = "disk"
    _con_args = ('center', '_norm_vec', '_radius', '_height')
    def __init__(self, center, normal, radius, height, fields=None,
                 pf=None, **kwargs):
        """
        By providing a *center*, a *normal*, a *radius* and a *height* we
        can define a cylinder of any proportion.  Only cells whose centers are
        within the cylinder will be selected.
        """
        AMR3DData.__init__(self, na.array(center), fields, pf, **kwargs)
        self._norm_vec = na.array(normal)/na.sqrt(na.dot(normal,normal))
        self.set_field_parameter("height_vector", self._norm_vec)
        self._height = height
        self._radius = radius
        self._d = -1.0 * na.dot(self._norm_vec, self.center)
        self._refresh_data()

    def _get_list_of_grids(self):
        H = na.sum(self._norm_vec.reshape((1,3,1)) * self.pf.h.gridCorners,
                   axis=1) + self._d
        D = na.sqrt(na.sum((self.pf.h.gridCorners -
                           self.center.reshape((1,3,1)))**2.0,axis=1))
        R = na.sqrt(D**2.0-H**2.0)
        self._grids = self.hierarchy.grids[
            ( (na.any(na.abs(H)<self._height,axis=0))
            & (na.any(R<self._radius,axis=0)
            & (na.logical_not((na.all(H>0,axis=0) | (na.all(H<0, axis=0)))) )
            ) ) ]
        self._grids = self.hierarchy.grids

    def _is_fully_enclosed(self, grid):
        corners = grid._corners.reshape((8,3,1))
        H = na.sum(self._norm_vec.reshape((1,3,1)) * corners,
                   axis=1) + self._d
        D = na.sqrt(na.sum((corners -
                           self.center.reshape((1,3,1)))**2.0,axis=1))
        R = na.sqrt(D**2.0-H**2.0)
        return (na.all(na.abs(H) < self._height, axis=0) \
            and na.all(R < self._radius, axis=0))

    @cache_mask
    def _get_cut_mask(self, grid):
        if self._is_fully_enclosed(grid):
            return True
        else:
            h = grid['x'] * self._norm_vec[0] \
              + grid['y'] * self._norm_vec[1] \
              + grid['z'] * self._norm_vec[2] \
              + self._d
            d = na.sqrt(
                (grid['x'] - self.center[0])**2.0
              + (grid['y'] - self.center[1])**2.0
              + (grid['z'] - self.center[2])**2.0
                )
            r = na.sqrt(d**2.0-h**2.0)
            cm = ( (na.abs(h) < self._height)
                 & (r < self._radius))
        return cm

class AMRRegionBase(AMR3DData):
    """
    AMRRegions are rectangular prisms of data.
    """
    _type_name = "region"
    _con_args = ('center', 'left_edge', 'right_edge')
    _dx_pad = 0.5
    def __init__(self, center, left_edge, right_edge, fields = None,
                 pf = None, **kwargs):
        """
        We create an object with a set of three *left_edge* coordinates,
        three *right_edge* coordinates, and a *center* that need not be the
        center.
        """
        AMR3DData.__init__(self, center, fields, pf, **kwargs)
        self.left_edge = left_edge
        self.right_edge = right_edge
        self._refresh_data()

    def _get_list_of_grids(self):
        self._grids, ind = self.pf.hierarchy.get_box_grids(self.left_edge,
                                                           self.right_edge)

    def _is_fully_enclosed(self, grid):
        return na.all( (grid._corners < self.right_edge)
                     & (grid._corners > self.left_edge))

    @cache_mask
    def _get_cut_mask(self, grid):
        if self._is_fully_enclosed(grid):
            return True
        else:
            dxp, dyp, dzp = self._dx_pad * grid.dds
            cm = ( (grid['x'] - dxp < self.right_edge[0])
                 & (grid['x'] + dxp > self.left_edge[0])
                 & (grid['y'] - dyp < self.right_edge[1])
                 & (grid['y'] + dyp > self.left_edge[1])
                 & (grid['z'] - dzp < self.right_edge[2])
                 & (grid['z'] + dzp > self.left_edge[2]) )
        return cm

class AMRRegionStrictBase(AMRRegionBase):
    """
    AMRRegion without any dx padding for cell selection
    """
    _dx_pad = 0.0

class AMRPeriodicRegionBase(AMR3DData):
    """
    AMRRegions are rectangular prisms of data.
    """
    _type_name = "periodic_region"
    _con_args = ('center', 'left_edge', 'right_edge')
    _dx_pad = 0.5
    def __init__(self, center, left_edge, right_edge, fields = None,
                 pf = None, **kwargs):
        """
        We create an object with a set of three *left_edge* coordinates,
        three *right_edge* coordinates, and a *center* that need not be the
        center.
        """
        AMR3DData.__init__(self, center, fields, pf, **kwargs)
        self.left_edge = na.array(left_edge)
        self.right_edge = na.array(right_edge)
        self._refresh_data()
        self.offsets = (na.mgrid[-1:1:3j,-1:1:3j,-1:1:3j] * \
                        (self.pf["DomainRightEdge"] -
                         self.pf["DomainLeftEdge"])[:,None,None,None])\
                       .transpose().reshape(27,3) # cached and in order

    def _get_list_of_grids(self):
        self._grids, ind = self.pf.hierarchy.get_periodic_box_grids(self.left_edge,
                                                                    self.right_edge)

    def _is_fully_enclosed(self, grid):
        for off_x, off_y, off_z in self.offsets:
            region_left = [self.left_edge[0]+off_x,
                           self.left_edge[1]+off_y,self.left_edge[2]+off_z]
            region_right = [self.right_edge[0]+off_x,
                            self.right_edge[1]+off_y,self.right_edge[2]+off_z]
            if (na.all((grid._corners <= region_right) &
                       (grid._corners >= region_left))):
                return True
        return False

    @cache_mask
    def _get_cut_mask(self, grid):
        if self._is_fully_enclosed(grid):
            return True
        else:
            cm = na.zeros(grid.ActiveDimensions,dtype='bool')
            dxp, dyp, dzp = self._dx_pad * grid.dds
            for off_x, off_y, off_z in self.offsets:
                cm = cm | ( (grid['x'] - dxp + off_x < self.right_edge[0])
                          & (grid['x'] + dxp + off_x > self.left_edge[0])
                          & (grid['y'] - dyp + off_y < self.right_edge[1])
                          & (grid['y'] + dyp + off_y > self.left_edge[1])
                          & (grid['z'] - dzp + off_z < self.right_edge[2])
                          & (grid['z'] + dzp + off_z > self.left_edge[2]) )
            return cm


class AMRPeriodicRegionStrictBase(AMRPeriodicRegionBase):
    """
    AMRPeriodicRegion without any dx padding for cell selection
    """
    _dx_pad = 0.0

class AMRGridCollection(AMR3DData):
    """
    An arbitrary selection of grids, within which we accept all points.
    """
    _type_name = "grid_collection"
    _con_args = ("center", "grid_list")
    def __init__(self, center, grid_list, fields = None,
                 pf = None, **kwargs):
        """
        By selecting an arbitrary *grid_list*, we can act on those grids.
        Child cells are not returned.
        """
        AMR3DData.__init__(self, center, fields, pf, **kwargs)
        self._grids = na.array(grid_list)

    def _get_list_of_grids(self):
        pass

    def _is_fully_enclosed(self, grid):
        return True

    @cache_mask
    def _get_cut_mask(self, grid):
        return na.ones(grid.ActiveDimensions, dtype='bool')

    def _get_point_indices(self, grid, use_child_mask=True):
        k = na.ones(grid.ActiveDimensions, dtype='bool')
        if use_child_mask:
            k[grid.child_indices] = False
        pointI = na.where(k == True)
        return pointI

class AMRSphereBase(AMR3DData):
    """
    A sphere of points
    """
    _type_name = "sphere"
    _con_args = ('center', 'radius')
    def __init__(self, center, radius, fields = None, pf = None, **kwargs):
        """
        The most famous of all the data objects, we define it via a
        *center* and a *radius*.
        """
        AMR3DData.__init__(self, center, fields, pf, **kwargs)
        if radius < self.hierarchy.get_smallest_dx():
            raise SyntaxError("Your radius is smaller than your finest cell!")
        self.set_field_parameter('radius',radius)
        self.radius = radius
        self.DW = self.pf["DomainRightEdge"] - self.pf["DomainLeftEdge"]
        self._refresh_data()

    def _get_list_of_grids(self, field = None):
        grids,ind = self.hierarchy.find_sphere_grids(self.center, self.radius)
        # Now we sort by level
        grids = grids.tolist()
        grids.sort(key=lambda x: (x.Level, x.LeftEdge[0], x.LeftEdge[1], x.LeftEdge[2]))
        self._grids = na.array(grids)

    def _is_fully_enclosed(self, grid):
        r = na.abs(grid._corners - self.center)
        r = na.minimum(r, na.abs(self.DW[None,:]-r))
        corner_radius = na.sqrt((r**2.0).sum(axis=1))
        return na.all(corner_radius <= self.radius)

    @restore_grid_state # Pains me not to decorate with cache_mask here
    def _get_cut_mask(self, grid, field=None):
        # We have the *property* center, which is not necessarily
        # the same as the field_parameter
        if self._is_fully_enclosed(grid):
            return True # We do not want child masking here
        if not isinstance(grid, (FakeGridForParticles, GridChildMaskWrapper)) \
           and grid.id in self._cut_masks:
            return self._cut_masks[grid.id]
        cm = ( (grid["RadiusCode"]<=self.radius) & grid.child_mask )
        if not isinstance(grid, (FakeGridForParticles, GridChildMaskWrapper)):
            self._cut_masks[grid.id] = cm
        return cm

class AMRCoveringGridBase(AMR3DData):
    """
    Covering grids represent fixed-resolution data over a given region.
    In order to achieve this goal -- for instance in order to obtain ghost
    zones -- grids up to and including the indicated level are included.
    No interpolation is done (as that would affect the 'power' on small
    scales) on the input data.
    """
    _spatial = True
    _type_name = "covering_grid"
    _con_args = ('level', 'left_edge', 'right_edge', 'ActiveDimensions')
    def __init__(self, level, left_edge, right_edge, dims, fields = None,
                 pf = None, num_ghost_zones = 0, use_pbar = True, **kwargs):
        """
        The data object returned will consider grids up to *level* in
        generating fixed resolution data between *left_edge* and *right_edge*
        that is *dims* (3-values) on a side.
        """
        AMR3DData.__init__(self, center=None, fields=fields, pf=pf, **kwargs)
        self.left_edge = na.array(left_edge)
        self.right_edge = na.array(right_edge)
        self.level = level
        self.ActiveDimensions = na.array(dims)
        dds = (self.right_edge-self.left_edge) \
              / self.ActiveDimensions
        self.dds = dds
        self.data["dx"] = dds[0]
        self.data["dy"] = dds[1]
        self.data["dz"] = dds[2]
        self._num_ghost_zones = num_ghost_zones
        self._use_pbar = use_pbar
        self._refresh_data()

    def _get_list_of_grids(self):
        if self._grids is not None: return
        if na.any(self.left_edge < self.pf["DomainLeftEdge"]) or \
           na.any(self.right_edge > self.pf["DomainRightEdge"]):
            grids,ind = self.pf.hierarchy.get_periodic_box_grids(
                            self.left_edge, self.right_edge)
            ind = slice(None)
        else:
            grids,ind = self.pf.hierarchy.get_box_grids(
                            self.left_edge, self.right_edge)
        level_ind = na.where(self.pf.hierarchy.gridLevels.ravel()[ind] <= self.level)
        sort_ind = na.argsort(self.pf.h.gridLevels.ravel()[ind][level_ind])
        self._grids = self.pf.hierarchy.grids[ind][level_ind][(sort_ind,)][::-1]

    def extract_region(self, indices):
        mylog.error("Sorry, dude, do it yourself, it's already in 3-D.")

    def _refresh_data(self):
        AMR3DData._refresh_data(self)
        self['dx'] = self.dds[0] * na.ones(self.ActiveDimensions, dtype='float64')
        self['dy'] = self.dds[1] * na.ones(self.ActiveDimensions, dtype='float64')
        self['dz'] = self.dds[2] * na.ones(self.ActiveDimensions, dtype='float64')

    def get_data(self, fields=None):
        if self._grids is None:
            self._get_list_of_grids()
        if fields is None:
            fields = self.fields[:]
        else:
            fields = ensure_list(fields)
        obtain_fields = []
        for field in fields:
            if self.data.has_key(field): continue
            if field not in self.hierarchy.field_list:
                try:
                    #print "Generating", field
                    self._generate_field(field)
                    continue
                except NeedsOriginalGrid, ngt_exception:
                    pass
            obtain_fields.append(field)
            self[field] = na.zeros(self.ActiveDimensions, dtype='float64') -999
        if len(obtain_fields) == 0: return
        mylog.debug("Getting fields %s from %s possible grids",
                   obtain_fields, len(self._grids))
        if self._use_pbar: pbar = \
                get_pbar('Searching grids for values ', len(self._grids))
        for i, grid in enumerate(self._grids):
            if self._use_pbar: pbar.update(i)
            self._get_data_from_grid(grid, obtain_fields)
            if not na.any(self[obtain_fields[0]] == -999): break
        if self._use_pbar: pbar.finish()
        if na.any(self[obtain_fields[0]] == -999):
            # and self.dx < self.hierarchy.grids[0].dx:
            print "COVERING PROBLEM", na.where(self[obtain_fields[0]]==-999)[0].size
            print na.where(self[obtain_fields[0]]==-999)
            raise KeyError
            
    def _generate_field(self, field):
        if self.pf.field_info.has_key(field):
            # First we check the validator; this might even raise!
            self.pf.field_info[field].check_available(self)
            self[field] = self.pf.field_info[field](self)
        else: # Can't find the field, try as it might
            raise exceptions.KeyError(field)

    def flush_data(self, field=None):
        """
        Any modifications made to the data in this object are pushed back
        to the originating grids, except the cells where those grids are both
        below the current level `and` have child cells.
        """
        self._get_list_of_grids()
        # We don't generate coordinates here.
        if field == None:
            fields_to_get = self.fields
        else:
            fields_to_get = ensure_list(field)
        for grid in self._grids:
            self._flush_data_to_grid(grid, fields_to_get)

    @restore_grid_state
    def _get_data_from_grid(self, grid, fields):
        ll = int(grid.Level == self.level)
        g_dx = grid.dds.ravel()
        c_dx = self.dds.ravel()
        g_fields = [grid[field] for field in ensure_list(fields)]
        c_fields = [self[field] for field in ensure_list(fields)]
        PointCombine.DataCubeRefine(
            grid.LeftEdge, g_dx, g_fields, grid.child_mask,
            self.left_edge, self.right_edge, c_dx, c_fields,
            ll, self.pf["DomainLeftEdge"], self.pf["DomainRightEdge"])

    def _flush_data_to_grid(self, grid, fields):
        ll = int(grid.Level == self.level)
        g_dx = grid.dds.ravel()
        c_dx = self.dds.ravel()
        g_fields = []
        for field in ensure_list(fields):
            if not grid.has_key(field): grid[field] = \
               na.zeros(grid.ActiveDimensions, dtype=self[field].dtype)
            g_fields.append(grid[field])
        c_fields = [self[field] for field in ensure_list(fields)]
        PointCombine.DataCubeReplace(
            grid.LeftEdge, g_dx, g_fields, grid.child_mask,
            self.left_edge, self.right_edge, c_dx, c_fields,
            ll, self.pf["DomainLeftEdge"], self.pf["DomainRightEdge"])

    @property
    def LeftEdge(self):
        return self.left_edge

    @property
    def RightEdge(self):
        return self.right_edge

class AMRSmoothedCoveringGridBase(AMRCoveringGridBase):
    _type_name = "smoothed_covering_grid"
    def __init__(self, *args, **kwargs):
        dlog2 = na.log10(kwargs['dims'])/na.log10(2)
        if not na.all(na.floor(dlog2) == na.ceil(dlog2)):
            mylog.warning("Must be power of two dimensions")
            #raise ValueError
        kwargs['num_ghost_zones'] = 0
        AMRCoveringGridBase.__init__(self, *args, **kwargs)
        if na.any(self.left_edge == self.pf["DomainLeftEdge"]):
            self.left_edge += self.dds
            self.ActiveDimensions -= 1
        if na.any(self.right_edge == self.pf["DomainRightEdge"]):
            self.right_edge -= self.dds
            self.ActiveDimensions -= 1

    def _get_list_of_grids(self):
        if na.any(self.left_edge - self.dds < self.pf["DomainLeftEdge"]) or \
           na.any(self.right_edge + self.dds > self.pf["DomainRightEdge"]):
            grids,ind = self.pf.hierarchy.get_periodic_box_grids(
                            self.left_edge - self.dds,
                            self.right_edge + self.dds)
            ind = slice(None)
        else:
            grids,ind = self.pf.hierarchy.get_box_grids(
                            self.left_edge - self.dds,
                            self.right_edge + self.dds)
        level_ind = na.where(self.pf.hierarchy.gridLevels.ravel()[ind] <= self.level)
        sort_ind = na.argsort(self.pf.h.gridLevels.ravel()[ind][level_ind])
        self._grids = self.pf.hierarchy.grids[ind][level_ind][(sort_ind,)]

    def _get_level_array(self, level, fields):
        fields = ensure_list(fields)
        # We assume refinement by a factor of two
        rf = float(self.pf["RefineBy"]**(self.level - level))
        dims = na.maximum(1,self.ActiveDimensions/rf) + 2
        dx = (self.right_edge-self.left_edge)/(dims-2)
        x,y,z = (na.mgrid[0:dims[0],0:dims[1],0:dims[2]].astype('float64')+0.5)\
              * dx[0]
        x += self.left_edge[0] - dx[0]
        y += self.left_edge[1] - dx[1]
        z += self.left_edge[2] - dx[2]
        offsets = [self['cd%s' % ax]*0.5 for ax in 'xyz']
        bounds = [self.left_edge[0]-offsets[0], self.right_edge[0]+offsets[0],
                  self.left_edge[1]-offsets[1], self.right_edge[1]+offsets[1],
                  self.left_edge[2]-offsets[2], self.right_edge[2]+offsets[2]]
        fake_grid = {'x':x,'y':y,'z':z,'dx':dx[0],'dy':dx[1],'dz':dx[2]}
        for ax in 'xyz': self['cd%s'%ax] = fake_grid['d%s'%ax]
        for field in fields:
            # Generate the new grid field
            if field in self.pf.field_info and self.pf.field_info[field].take_log:
                interpolator = TrilinearFieldInterpolator(
                                na.log10(self[field]), bounds, ['x','y','z'],
                                truncate = True)
                self[field] = 10**interpolator(fake_grid)
            else:
                interpolator = TrilinearFieldInterpolator(
                                self[field], bounds, ['x','y','z'],
                                truncate = True)
                self[field] = interpolator(fake_grid)
        return fake_grid

    def get_data(self, field=None):
        self._get_list_of_grids()
        # We don't generate coordinates here.
        if field == None:
            fields_to_get = self.fields
        else:
            fields_to_get = ensure_list(field)
        for field in fields_to_get:
            grid_count = 0
            if self.data.has_key(field):
                continue
            mylog.debug("Getting field %s from %s possible grids",
                       field, len(self._grids))
            self[field] = na.zeros(self.ActiveDimensions, dtype='float64') -999
            if self._use_pbar: pbar = \
                    get_pbar('Searching grids for values ', len(self._grids))
            # How do we find out the root grid base dx?
            idims = na.array([3,3,3])
            dx = na.minimum((self.right_edge-self.left_edge)/(idims-2),
                            self.pf.h.grids[0]['dx'])
            idims = na.floor((self.right_edge-self.left_edge)/dx) + 2
            for ax in 'xyz': self['cd%s'%ax] = dx[0]
            self[field] = na.zeros(idims,dtype='float64')-999
            for level in range(self.level+1):
                for grid in self.select_grids(level):
                    if self._use_pbar: pbar.update(grid_count)
                    self._get_data_from_grid(grid, field)
                    grid_count += 1
                if level < self.level: self._get_level_array(level+1, field)
            self[field] = self[field][1:-1,1:-1,1:-1]
            if self._use_pbar: pbar.finish()
        
    @restore_grid_state
    def _get_data_from_grid(self, grid, fields):
        fields = ensure_list(fields)
        g_dx = grid.dds
        c_dx = na.array([self['cdx'],self['cdy'],self['cdz']])
        g_fields = [grid[field] for field in fields]
        c_fields = [self[field] for field in fields]
        total = PointCombine.DataCubeRefine(
            grid.LeftEdge, g_dx, g_fields, grid.child_mask,
            self.left_edge-c_dx, self.right_edge+c_dx,
            c_dx, c_fields,
            1, self.pf["DomainLeftEdge"], self.pf["DomainRightEdge"])

    def flush_data(self, *args, **kwargs):
        raise KeyError("Can't do this")


class EnzoOrthoRayBase(AMROrthoRayBase): pass
class EnzoRayBase(AMRRayBase): pass
class EnzoSliceBase(AMRSliceBase): pass
class EnzoCuttingPlaneBase(AMRCuttingPlaneBase): pass
class EnzoProjBase(AMRProjBase): pass
class EnzoCylinderBase(AMRCylinderBase): pass
class EnzoRegionBase(AMRRegionBase): pass
class EnzoPeriodicRegionBase(AMRPeriodicRegionBase): pass
class EnzoGridCollection(AMRGridCollection): pass
class EnzoSphereBase(AMRSphereBase): pass
class EnzoCoveringGrid(AMRCoveringGridBase): pass
class EnzoSmoothedCoveringGrid(AMRSmoothedCoveringGridBase): pass

def _reconstruct_object(*args, **kwargs):
    pfid = args[0]
    dtype = args[1]
    field_parameters = args[-1]
    # will be much nicer when we can do pfid, *a, fp = args
    args, new_args = args[2:-1], []
    for arg in args:
        if iterable(arg) and len(arg) == 2 \
           and not isinstance(arg, types.DictType) \
           and isinstance(arg[1], AMRData):
            new_args.append(arg[1])
        else: new_args.append(arg)
    pfs = ParameterFileStore()
    pf = pfs.get_pf_hash(pfid)
    cls = getattr(pf.h, dtype)
    obj = cls(*new_args)
    obj.field_parameters.update(field_parameters)
    return pf, obj
