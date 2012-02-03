"""
Data containers based on geometric selection

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Author: Britton Smith <Britton.Smith@colorado.edu>
Affiliation: University of Colorado at Boulder
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

import types
import numpy as na
from exceptions import ValueError, SyntaxError

from yt.funcs import *
from yt.utilities.amr_utils import \
    VoxelTraversal, planar_points_in_volume, find_grids_in_inclined_box, \
    grid_points_in_volume
from .data_containers import \
    YTSelectionContainer1D, YTSelectionContainer2D, YTSelectionContainer3D, \
    restore_grid_state, cache_mask, cache_point_indices, cache_vc_data, \
    FakeGridForParticles
from yt.geometry.selection_routines import \
    ortho_ray_grids, ray_grids, slice_grids, cutting_plane_grids
from yt.data_objects.derived_quantities import \
    DerivedQuantityCollection, GridChildMaskWrapper
from yt.utilities.definitions import \
    x_dict, y_dict, axis_names
from yt.utilities.exceptions import YTSphereTooSmall
from yt.utilities.linear_interpolators import TrilinearFieldInterpolator

class YTOrthoRayBase(YTSelectionContainer1D):
    _key_fields = ['x','y','z','dx','dy','dz']
    _type_name = "ortho_ray"
    _con_args = ('axis', 'coords')
    def __init__(self, axis, coords, fields=None, pf=None, **kwargs):
        """
        This is an orthogonal ray cast through the entire domain, at a specific
        coordinate.

        This object is typically accessed through the `ortho_ray` object that
        hangs off of hierarchy objects.  The resulting arrays have their
        dimensionality reduced to one, and an ordered list of points at an
        (x,y) tuple along `axis` are available.

        Parameters
        ----------
        axis : int
            The axis along which to cast the ray.  Can be 0, 1, or 2 for x, y, z.
        coords : tuple of floats
            The (plane_x, plane_y) coordinates at which to cast the ray.  Note
            that this is in the plane coordinates: so if you are casting along
            x, this will be (y,z).  If you are casting along y, this will be
            (x,z).  If you are casting along z, this will be (x,y).
        fields : list of strings, optional
            If you want the object to pre-retrieve a set of fields, supply them
            here.  This is not necessary.
        kwargs : dict of items
            Any additional values are passed as field parameters that can be
            accessed by generated fields.

        Examples
        --------

        >>> pf = load("RedshiftOutput0005")
        >>> oray = pf.h.ortho_ray(0, (0.2, 0.74))
        >>> print oray["Density"]
        """
        YTSelectionContainer1D.__init__(self, pf, fields, **kwargs)
        self.axis = axis
        self.px_ax = x_dict[self.axis]
        self.py_ax = y_dict[self.axis]
        self.px_dx = 'd%s'%(axis_names[self.px_ax])
        self.py_dx = 'd%s'%(axis_names[self.py_ax])
        self.px, self.py = coords
        self.sort_by = axis_names[self.axis]
        self._refresh_data()

    @property
    def coords(self):
        return (self.px, self.py)

    def _get_list_of_grids(self):
        gi = ortho_ray_grids(self,
                self.hierarchy.grid_left_edge,
                self.hierarchy.grid_right_edge)
        self._grids = self.hierarchy.grids[gi]

    @restore_grid_state
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


class YTRayBase(YTSelectionContainer1D):
    _type_name = "ray"
    _con_args = ('start_point', 'end_point')
    sort_by = 't'
    def __init__(self, start_point, end_point, fields=None, pf=None, **kwargs):
        """
        This is an arbitrarily-aligned ray cast through the entire domain, at a
        specific coordinate.

        This object is typically accessed through the `ray` object that hangs
        off of hierarchy objects.  The resulting arrays have their
        dimensionality reduced to one, and an ordered list of points at an
        (x,y) tuple along `axis` are available, as is the `t` field, which
        corresponds to a unitless measurement along the ray from start to
        end.

        Parameters
        ----------
        start_point : array-like set of 3 floats
            The place where the ray starts.
        end_point : array-like set of 3 floats
            The place where the ray ends.
        fields : list of strings, optional
            If you want the object to pre-retrieve a set of fields, supply them
            here.  This is not necessary.
        kwargs : dict of items
            Any additional values are passed as field parameters that can be
            accessed by generated fields.

        Examples
        --------

        >>> pf = load("RedshiftOutput0005")
        >>> ray = pf.h._ray((0.2, 0.74, 0.11), (0.4, 0.91, 0.31))
        >>> print ray["Density"], ray["t"], ray["dts"]
        """
        YTSelectionContainer1D.__init__(self, pf, fields, **kwargs)
        self.start_point = na.array(start_point, dtype='float64')
        self.end_point = na.array(end_point, dtype='float64')
        self.vec = self.end_point - self.start_point
        #self.vec /= na.sqrt(na.dot(self.vec, self.vec))
        self._set_center(self.start_point)
        self.set_field_parameter('center', self.start_point)
        self._dts, self._ts = {}, {}
        #self._refresh_data()

    def _get_list_of_grids(self):
        gi = ray_grids(self,
                self.hierarchy.grid_left_edge,
                self.hierarchy.grid_right_edge)
        self._grids = self.hierarchy.grids[gi]

    @restore_grid_state
    def _get_data_from_grid(self, grid, field):
        mask = na.logical_and(self._get_cut_mask(grid),
                              grid.child_mask)
        if field == 'dts': return self._dts[grid.id][mask]
        if field == 't': return self._ts[grid.id][mask]
        gf = grid[field]
        if not iterable(gf):
            gf = gf * na.ones(grid.child_mask.shape)
        return gf[mask]

    @cache_mask
    def _get_cut_mask(self, grid):
        mask = na.zeros(grid.ActiveDimensions, dtype='int')
        dts = na.zeros(grid.ActiveDimensions, dtype='float64')
        ts = na.zeros(grid.ActiveDimensions, dtype='float64')
        VoxelTraversal(mask, ts, dts, grid.LeftEdge, grid.RightEdge,
                       grid.dds, self.center, self.vec)
        self._dts[grid.id] = na.abs(dts)
        self._ts[grid.id] = na.abs(ts)
        return mask


class YTSliceBase(YTSelectionContainer2D):
    _top_node = "/Slices"
    _type_name = "slice"
    _con_args = ('axis', 'coord')
    #@time_execution
    def __init__(self, axis, coord, fields = None, center=None, pf=None,
                 node_name = False, **kwargs):
        """
        This is a data object corresponding to a slice through the simulation
        domain.

        This object is typically accessed through the `slice` object that hangs
        off of hierarchy objects.  AMRSlice is an orthogonal slice through the
        data, taking all the points at the finest resolution available and then
        indexing them.  It is more appropriately thought of as a slice
        'operator' than an object, however, as its field and coordinate can
        both change.

        Parameters
        ----------
        axis : int
            The axis along which to slice.  Can be 0, 1, or 2 for x, y, z.
        coord : float
            The coordinate along the axis at which to slice.  This is in
            "domain" coordinates.
        fields : list of strings, optional
            If you want the object to pre-retrieve a set of fields, supply them
            here.  This is not necessary.
        center : array_like, optional
            The 'center' supplied to fields that use it.  Note that this does
            not have to have `coord` as one value.  Strictly optional.
        node_name: string, optional
            The node in the .yt file to find or store this slice at.  Should
            probably not be used.
        kwargs : dict of items
            Any additional values are passed as field parameters that can be
            accessed by generated fields.

        Examples
        --------

        >>> pf = load("RedshiftOutput0005")
        >>> slice = pf.h.slice(0, 0.25)
        >>> print slice["Density"]
        """
        YTSelectionContainer2D.__init__(self, axis, fields, pf, **kwargs)
        self._set_center(center)
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
            dx = self.hierarchy.select_grids(level)[0].dds[self.axis]
            self.coord += dx * val
        else:
            raise ValueError(val)
        self._refresh_data()

    def _generate_coords(self):
        points = []
        for grid in self._get_grids():
            points.append(self._generate_grid_coords(grid))
        if len(points) == 0:
            points = None
            t = self.comm.par_combine_object(None, datatype="array", op="cat")
        else:
            points = na.concatenate(points)
            # We have to transpose here so that _par_combine_object works
            # properly, as it and the alltoall assume the long axis is the last
            # one.
            t = self.comm.par_combine_object(points.transpose(),
                        datatype="array", op="cat")
        self['px'] = t[0,:]
        self['py'] = t[1,:]
        self['pz'] = t[2,:]
        self['pdx'] = t[3,:]
        self['pdy'] = t[4,:]
        self['pdz'] = t[3,:] # Does not matter!

        # Now we set the *actual* coordinates
        self[axis_names[x_dict[self.axis]]] = t[0,:]
        self[axis_names[y_dict[self.axis]]] = t[1,:]
        self[axis_names[self.axis]] = t[2,:]

        self.ActiveDimensions = (t.shape[1], 1, 1)

    def _get_list_of_grids(self):
        gi = slice_grids(self,
                self.hierarchy.grid_left_edge,
                self.hierarchy.grid_right_edge)
        self._grids = self.hierarchy.grids[gi]

    def __cut_mask_child_mask(self, grid):
        mask = grid.child_mask.copy()
        return mask

    def _generate_grid_coords(self, grid):
        xaxis = x_dict[self.axis]
        yaxis = y_dict[self.axis]
        ds, dx, dy = grid.dds[self.axis], grid.dds[xaxis], grid.dds[yaxis]
        sl_ind = int((self.coord-self.pf.domain_left_edge[self.axis])/ds) - \
                     grid.get_global_startindex()[self.axis]
        sl = [slice(None), slice(None), slice(None)]
        sl[self.axis] = slice(sl_ind, sl_ind + 1)
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
        sl_ind = int((self.coord-self.pf.domain_left_edge[self.axis])/dx) - \
                     grid.get_global_startindex()[self.axis]
        sl = [slice(None), slice(None), slice(None)]
        sl[self.axis] = slice(sl_ind, sl_ind + 1)
        sl = tuple(sl)
        if self.pf.field_info.has_key(field) and self.pf.field_info[field].particle_type:
            return grid[field]
        elif field in self.pf.field_info and self.pf.field_info[field].not_in_all:
            dv = grid[field][sl]
        elif not grid.has_key(field):
            conv_factor = 1.0
            if self.pf.field_info.has_key(field):
                conv_factor = self.pf.field_info[field]._convert_function(self)
            dv = self.hierarchy.io._read_data_slice(grid, field, self.axis, sl_ind) * conv_factor
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

    def __get_quantities(self):
        if self.__quantities is None:
            self.__quantities = DerivedQuantityCollection(self)
        return self.__quantities
    __quantities = None
    quantities = property(__get_quantities)


class YTCuttingPlaneBase(YTSelectionContainer2D):
    _plane = None
    _top_node = "/CuttingPlanes"
    _key_fields = YTSelectionContainer2D._key_fields + ['pz','pdz']
    _type_name = "cutting"
    _con_args = ('normal', 'center')
    def __init__(self, normal, center, fields = None, node_name = None,
                 **kwargs):
        """
        This is a data object corresponding to an oblique slice through the
        simulation domain.

        This object is typically accessed through the `cutting` object
        that hangs off of hierarchy objects.  AMRCuttingPlane is an oblique
        plane through the data, defined by a normal vector and a coordinate.
        It attempts to guess an 'up' vector, which cannot be overridden, and
        then it pixelizes the appropriate data onto the plane without
        interpolation.

        Parameters
        ----------
        normal : array_like
            The vector that defines the desired plane.  For instance, the
            angular momentum of a sphere.
        center : array_like, optional
            The center of the cutting plane.
        fields : list of strings, optional
            If you want the object to pre-retrieve a set of fields, supply them
            here.  This is not necessary.
        node_name: string, optional
            The node in the .yt file to find or store this slice at.  Should
            probably not be used.
        kwargs : dict of items
            Any additional values are passed as field parameters that can be
            accessed by generated fields.

        Notes
        -----

        This data object in particular can be somewhat expensive to create.
        It's also important to note that unlike the other 2D data objects, this
        oject provides px, py, pz, as some cells may have a height from the
        plane.

        Examples
        --------

        >>> pf = load("RedshiftOutput0005")
        >>> cp = pf.h.cutting([0.1, 0.2, -0.9], [0.5, 0.42, 0.6])
        >>> print cp["Density"]
        """
        YTSelectionContainer2D.__init__(self, 4, fields, **kwargs)
        self._set_center(center)
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
        gridi = cutting_plane_grids(self, self.pf.h.grid_left_edge,
                                          self.pf.h.grid_right_edge)
        self._grids = self.hierarchy.grids[gridi.astype("bool")]

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
        if len(points) == 0: points = None
        else: points = na.concatenate(points)
        t = self.comm.par_combine_object(points, datatype="array", op="cat")
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
        cen_name = ("%s" % (self.center,)).replace(" ","_")[1:-1]
        L_name = ("%s" % self._norm_vec).replace(" ","_")[1:-1]
        return "%s/c%s_L%s" % \
            (self._top_node, cen_name, L_name)

    def to_frb(self, width, resolution):
        r"""This function returns an ObliqueFixedResolutionBuffer generated
        from this object.

        An ObliqueFixedResolutionBuffer is an object that accepts a
        variable-resolution 2D object and transforms it into an NxM bitmap that
        can be plotted, examined or processed.  This is a convenience function
        to return an FRB directly from an existing 2D data object.  Unlike the
        corresponding to_frb function for other YTSelectionContainer2D objects, this does
        not accept a 'center' parameter as it is assumed to be centered at the
        center of the cutting plane.

        Parameters
        ----------
        width : width specifier
            This can either be a floating point value, in the native domain
            units of the simulation, or a tuple of the (value, unit) style.
            This will be the width of the FRB.
        resolution : int or tuple of ints
            The number of pixels on a side of the final FRB.

        Returns
        -------
        frb : :class:`~yt.visualization.fixed_resolution.ObliqueFixedResolutionBuffer`
            A fixed resolution buffer, which can be queried for fields.

        Examples
        --------

        >>> v, c = pf.h.find_max("Density")
        >>> sp = pf.h.sphere(c, (100.0, 'au'))
        >>> L = sp.quantities["AngularMomentumVector"]()
        >>> cutting = pf.h.cutting(L, c)
        >>> frb = cutting.to_frb( (1.0, 'pc'), 1024)
        >>> write_image(na.log10(frb["Density"]), 'density_1pc.png')
        """
        if iterable(width):
            w, u = width
            width = w/self.pf[u]
        if not iterable(resolution):
            resolution = (resolution, resolution)
        from yt.visualization.fixed_resolution import ObliqueFixedResolutionBuffer
        bounds = (-width/2.0, width/2.0, -width/2.0, width/2.0)
        frb = ObliqueFixedResolutionBuffer(self, bounds, resolution)
        return frb


class YTFixedResCuttingPlaneBase(YTSelectionContainer2D):
    """
    YTFixedResCuttingPlaneBase is an oblique plane through the data,
    defined by a normal vector and a coordinate.  It trilinearly
    interpolates the data to a fixed resolution slice.  It differs from
    the other data objects as it doesn't save the grid data, only the
    interpolated data.
    """
    _top_node = "/FixedResCuttingPlanes"
    _type_name = "fixed_res_cutting"
    _con_args = ('normal', 'center', 'width', 'dims')
    def __init__(self, normal, center, width, dims, fields = None,
                 node_name = None, **kwargs):
        """
        The fixed resolution Cutting Plane slices at an oblique angle,
        where we use the *normal* vector at the *center* to define the
        viewing plane.  The plane is *width* units wide.  The 'up'
        direction is guessed at automatically if not given.
        """
        #
        # Taken from Cutting Plane
        #
        YTSelectionContainer2D.__init__(self, 4, fields, **kwargs)
        self._set_center(center)
        self.width = width
        self.dims = dims
        self.dds = self.width / self.dims
        self.bounds = na.array([0.0,1.0,0.0,1.0])

        self.set_field_parameter('center', center)
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

        # Calculate coordinates of each pixel
        _co = self.dds * \
              (na.mgrid[-self.dims/2 : self.dims/2,
                        -self.dims/2 : self.dims/2] + 0.5)
        self._coord = self.center + na.outer(_co[0,:,:], self._x_vec) + \
                      na.outer(_co[1,:,:], self._y_vec)
        self._pixelmask = na.ones(self.dims*self.dims, dtype='int8')

        if node_name is False:
            self._refresh_data()
        else:
            if node_name is True: self._deserialize()
            else: self._deserialize(node_name)

    @property
    def normal(self):
        return self._norm_vec

    def _get_list_of_grids(self):
        # Just like the Cutting Plane but restrict the grids to be
        # within width/2 of the center.
        vertices = self.hierarchy.gridCorners
        # Shape = (8,3,n_grid)
        D = na.sum(self._norm_vec.reshape((1,3,1)) * vertices, axis=1) + self._d
        valid_grids = na.where(na.logical_not(na.all(D<0,axis=0) |
                                              na.all(D>0,axis=0) ))[0]
        # Now restrict these grids to a rect. prism that bounds the slice
        sliceCorners = na.array([ \
            self.center + 0.5*self.width * (+self._x_vec + self._y_vec),
            self.center + 0.5*self.width * (+self._x_vec - self._y_vec),
            self.center + 0.5*self.width * (-self._x_vec - self._y_vec),
            self.center + 0.5*self.width * (-self._x_vec + self._y_vec) ])
        sliceLeftEdge = sliceCorners.min(axis=0)
        sliceRightEdge = sliceCorners.max(axis=0)
        # Check for bounding box and grid overlap
        leftOverlap = na.less(self.hierarchy.gridLeftEdge[valid_grids],
                              sliceRightEdge).all(axis=1)
        rightOverlap = na.greater(self.hierarchy.gridRightEdge[valid_grids],
                                  sliceLeftEdge).all(axis=1)
        self._grids = self.hierarchy.grids[valid_grids[
            na.where(leftOverlap & rightOverlap)]]
        self._grids = self._grids[::-1]

    def _generate_coords(self):
        self['px'] = self._coord[:,0].ravel()
        self['py'] = self._coord[:,1].ravel()
        self['pz'] = self._coord[:,2].ravel()
        self['pdx'] = self.dds * 0.5
        self['pdy'] = self.dds * 0.5
        #self['pdz'] = self.dds * 0.5

    def _get_data_from_grid(self, grid, field):
        if not self.pf.field_info[field].particle_type:
            pointI = self._get_point_indices(grid)
            if len(pointI) == 0: return
            vc = self._calc_vertex_centered_data(grid, field)
            bds = na.array(zip(grid.LeftEdge,
                               grid.RightEdge)).ravel()
            interp = TrilinearFieldInterpolator(vc, bds, ['x', 'y', 'z'])
            self[field][pointI] = interp( \
                dict(x=self._coord[pointI,0],
                     y=self._coord[pointI,1],
                     z=self._coord[pointI,2])).ravel()

            # Mark these pixels to speed things up
            self._pixelmask[pointI] = 0

            return
        else:
            raise SyntaxError("Making a fixed resolution slice with "
                              "particles isn't supported yet.")

    def reslice(self, normal, center, width):

        # Cleanup
        del self._coord
        del self._pixelmask

        self.center = center
        self.width = width
        self.dds = self.width / self.dims
        self.set_field_parameter('center', center)
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
        self.set_field_parameter('cp_x_vec',self._x_vec)
        self.set_field_parameter('cp_y_vec',self._y_vec)
        self.set_field_parameter('cp_z_vec',self._norm_vec)
        # Calculate coordinates of each pixel
        _co = self.dds * \
              (na.mgrid[-self.dims/2 : self.dims/2,
                        -self.dims/2 : self.dims/2] + 0.5)

        self._coord = self.center + na.outer(_co[0,:,:], self._x_vec) + \
                      na.outer(_co[1,:,:], self._y_vec)
        self._pixelmask = na.ones(self.dims*self.dims, dtype='int8')

        self._refresh_data()
        return

    #@time_execution
    def get_data(self, fields = None):
        """
        Iterates over the list of fields and generates/reads them all.
        """
        self._get_list_of_grids()
        if not self.has_key('pdx'):
            self._generate_coords()
        if fields == None:
            fields_to_get = self.fields[:]
        else:
            fields_to_get = ensure_list(fields)
        temp_data = {}
        _size = self.dims * self.dims
        for field in fields_to_get:
            if self.field_data.has_key(field): continue
            if field not in self.hierarchy.field_list:
                if self._generate_field(field):
                    continue # A "True" return means we did it
            if not self._vc_data.has_key(field):
                self._vc_data[field] = {}
            self[field] = na.zeros(_size, dtype='float64')
            for grid in self._get_grids():
                self._get_data_from_grid(grid, field)
            self[field] = self.comm.mpi_allreduce(\
                self[field], op='sum').reshape([self.dims]*2).transpose()

    def interpolate_discretize(self, *args, **kwargs):
        pass

    @cache_vc_data
    def _calc_vertex_centered_data(self, grid, field):
        #return grid.retrieve_ghost_zones(1, field, smoothed=False)
        return grid.get_vertex_centered_data(field)

    def _get_point_indices(self, grid):
        if self._pixelmask.max() == 0: return []
        k = planar_points_in_volume(self._coord, self._pixelmask,
                                    grid.LeftEdge, grid.RightEdge,
                                    grid.child_mask, just_one(grid['dx']))
        return k

    def _gen_node_name(self):
        cen_name = ("%s" % (self.center,)).replace(" ","_")[1:-1]
        L_name = ("%s" % self._norm_vec).replace(" ","_")[1:-1]
        return "%s/c%s_L%s" % \
            (self._top_node, cen_name, L_name)


class YTCylinderBase(YTSelectionContainer3D):
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
        YTSelectionContainer3D.__init__(self, center, fields, pf, **kwargs)
        self._norm_vec = na.array(normal)/na.sqrt(na.dot(normal,normal))
        self.set_field_parameter("height_vector", self._norm_vec)
        self._height = fix_length(height, self.pf)
        self._radius = fix_length(radius, self.pf)
        self._d = -1.0 * na.dot(self._norm_vec, self.center)
        self._refresh_data()

    def _get_list_of_grids(self):
        H = na.sum(self._norm_vec.reshape((1,3,1)) * self.pf.h.grid_corners,
                   axis=1) + self._d
        D = na.sqrt(na.sum((self.pf.h.grid_corners -
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
            cm = ( (na.abs(h) <= self._height)
                 & (r <= self._radius))
        return cm


class YTInclinedBoxBase(YTSelectionContainer3D):
    _type_name="inclined_box"
    _con_args = ('origin','box_vectors')

    def __init__(self, origin, box_vectors, fields=None,
                 pf=None, **kwargs):
        """
        A rectangular prism with arbitrary alignment to the computational
        domain.  *origin* is the origin of the box, while *box_vectors* is an
        array of ordering [ax, ijk] that describes the three vectors that
        describe the box.  No checks are done to ensure that the box satisfies
        a right-hand rule, but if it doesn't, behavior is undefined.
        """
        self.origin = na.array(origin)
        self.box_vectors = na.array(box_vectors, dtype='float64')
        self.box_lengths = (self.box_vectors**2.0).sum(axis=1)**0.5
        center = origin + 0.5*self.box_vectors.sum(axis=0)
        YTSelectionContainer3D.__init__(self, center, fields, pf, **kwargs)
        self._setup_rotation_parameters()
        self._refresh_data()

    def _setup_rotation_parameters(self):
        xv = self.box_vectors[0,:]
        yv = self.box_vectors[1,:]
        zv = self.box_vectors[2,:]
        self._x_vec = xv / na.sqrt(na.dot(xv, xv))
        self._y_vec = yv / na.sqrt(na.dot(yv, yv))
        self._z_vec = zv / na.sqrt(na.dot(zv, zv))
        self._rot_mat = na.array([self._x_vec,self._y_vec,self._z_vec])
        self._inv_mat = na.linalg.pinv(self._rot_mat)

    def _get_list_of_grids(self):
        if self._grids is not None: return
        GLE = self.pf.h.grid_left_edge
        GRE = self.pf.h.grid_right_edge
        goodI = find_grids_in_inclined_box(self.box_vectors, self.center,
                                           GLE, GRE)
        cgrids = self.pf.h.grids[goodI.astype('bool')]
       # find_grids_in_inclined_box seems to be broken.
        cgrids = self.pf.h.grids[:]
        grids = []
        for i,grid in enumerate(cgrids):
            v = grid_points_in_volume(self.box_lengths, self.origin,
                                      self._rot_mat, grid.LeftEdge,
                                      grid.RightEdge, grid.dds,
                                      grid.child_mask, 1)
            if v: grids.append(grid)
        self._grids = na.empty(len(grids), dtype='object')
        for gi, g in enumerate(grids): self._grids[gi] = g


    def _is_fully_enclosed(self, grid):
        # This should be written at some point.
        # We'd rotate all eight corners into the space of the box, then check to
        # see if all are enclosed.
        return False

    def _get_cut_mask(self, grid):
        if self._is_fully_enclosed(grid):
            return True
        pm = na.zeros(grid.ActiveDimensions, dtype='int32')
        grid_points_in_volume(self.box_lengths, self.origin,
                              self._rot_mat, grid.LeftEdge,
                              grid.RightEdge, grid.dds, pm, 0)
        return pm


class YTRegionBase(YTSelectionContainer3D):
    """
    AMRRegions are rectangular prisms of data.
    """
    _type_name = "region"
    _con_args = ('center', 'left_edge', 'right_edge')
    _dx_pad = 0.5
    def __init__(self, center, left_edge, right_edge, fields = None,
                 pf = None, **kwargs):
        """A 3D region of data with an arbitrary center.

        Takes an array of three *left_edge* coordinates, three
        *right_edge* coordinates, and a *center* that can be anywhere
        in the domain. If the selected region extends past the edges
        of the domain, no data will be found there, though the
        object's `left_edge` or `right_edge` are not modified.

        Parameters
        ----------
        center : array_like
            The center of the region
        left_edge : array_like
            The left edge of the region
        right_edge : array_like
            The right edge of the region
        """
        YTSelectionContainer3D.__init__(self, center, fields, pf, **kwargs)
        self.left_edge = left_edge
        self.right_edge = right_edge
        self._refresh_data()

class YTGridCollectionBase(YTSelectionContainer3D):
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
        YTSelectionContainer3D.__init__(self, center, fields, pf, **kwargs)
        self._grids = na.array(grid_list)
        self.grid_list = self._grids

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


class YTSphereBase(YTSelectionContainer3D):
    """
    A sphere of points
    """
    _type_name = "sphere"
    _con_args = ('center', 'radius')
    def __init__(self, center, radius, fields = None, pf = None, **kwargs):
        """A sphere f points defined by a *center* and a *radius*.

        Parameters
        ----------
        center : array_like
            The center of the sphere.
        radius : float
            The radius of the sphere.

        Examples
        --------
        >>> pf = load("DD0010/moving7_0010")
        >>> c = [0.5,0.5,0.5]
        >>> sphere = pf.h.sphere(c,1.*pf['kpc'])
        """
        YTSelectionContainer3D.__init__(self, center, fields, pf, **kwargs)
        # Unpack the radius, if necessary
        radius = fix_length(radius, self.pf)
        if radius < self.hierarchy.get_smallest_dx():
            raise YTSphereTooSmall(pf, radius, self.hierarchy.get_smallest_dx())
        self.set_field_parameter('radius',radius)
        self.radius = radius
        self.DW = self.pf.domain_right_edge - self.pf.domain_left_edge
        self._refresh_data()
