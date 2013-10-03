"""
Data containers based on geometric selection




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import types
import numpy as np
from exceptions import ValueError, SyntaxError

from yt.funcs import *
from yt.utilities.lib import \
    VoxelTraversal, planar_points_in_volume, find_grids_in_inclined_box, \
    grid_points_in_volume
from yt.utilities.lib.alt_ray_tracers import clyindrical_ray_trace
from yt.utilities.orientation import Orientation
from .data_containers import \
    YTSelectionContainer1D, YTSelectionContainer2D, YTSelectionContainer3D
from yt.data_objects.derived_quantities import \
    DerivedQuantityCollection
from yt.utilities.definitions import \
    x_dict, y_dict, axis_names
from yt.utilities.exceptions import YTSphereTooSmall
from yt.utilities.linear_interpolators import TrilinearFieldInterpolator
from yt.utilities.minimal_representation import \
    MinimalSliceData
from yt.utilities.math_utils import get_rotation_matrix

class YTOrthoRayBase(YTSelectionContainer1D):
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
    _key_fields = ['x','y','z','dx','dy','dz']
    _type_name = "ortho_ray"
    _con_args = ('axis', 'coords')
    def __init__(self, axis, coords, pf=None, field_parameters=None):
        super(YTOrthoRayBase, self).__init__(pf, field_parameters)
        self.axis = axis
        self.px_ax = x_dict[self.axis]
        self.py_ax = y_dict[self.axis]
        self.px_dx = 'd%s'%(axis_names[self.px_ax])
        self.py_dx = 'd%s'%(axis_names[self.py_ax])
        self.px, self.py = coords
        self.sort_by = axis_names[self.axis]

    @property
    def coords(self):
        return (self.px, self.py)

class YTRayBase(YTSelectionContainer1D):
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
    _type_name = "ray"
    _con_args = ('start_point', 'end_point')
    _container_fields = ("t", "dts")
    def __init__(self, start_point, end_point, pf=None, field_parameters=None):
        super(YTRayBase, self).__init__(pf, field_parameters)
        self.start_point = np.array(start_point, dtype='float64')
        self.end_point = np.array(end_point, dtype='float64')
        self.vec = self.end_point - self.start_point
        #self.vec /= np.sqrt(np.dot(self.vec, self.vec))
        self._set_center(self.start_point)
        self.set_field_parameter('center', self.start_point)
        self._dts, self._ts = None, None

    def _generate_container_field(self, field):
        if self._current_chunk is None:
            self.hierarchy._identify_base_chunk(self)
        if field == "dts":
            return self._current_chunk.dtcoords
        elif field == "t":
            return self._current_chunk.tcoords
        else:
            raise KeyError(field)

class YTSliceBase(YTSelectionContainer2D):
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
    _top_node = "/Slices"
    _type_name = "slice"
    _con_args = ('axis', 'coord')
    _container_fields = ("px", "py", "pdx", "pdy")

    def __init__(self, axis, coord, center=None, pf=None,
                 field_parameters = None):
        YTSelectionContainer2D.__init__(self, axis, pf, field_parameters)
        self._set_center(center)
        self.coord = coord

    def _generate_container_field(self, field):
        if self._current_chunk is None:
            self.hierarchy._identify_base_chunk(self)
        if field == "px":
            return self._current_chunk.fcoords[:,x_dict[self.axis]]
        elif field == "py":
            return self._current_chunk.fcoords[:,y_dict[self.axis]]
        elif field == "pdx":
            return self._current_chunk.fwidth[:,x_dict[self.axis]] * 0.5
        elif field == "pdy":
            return self._current_chunk.fwidth[:,y_dict[self.axis]] * 0.5
        else:
            raise KeyError(field)

    @property
    def _mrep(self):
        return MinimalSliceData(self)

    def hub_upload(self):
        self._mrep.upload()

    def to_pw(self, fields=None, center='c', width=None, axes_unit=None, 
               origin='center-window'):
        r"""Create a :class:`~yt.visualization.plot_window.PWViewerMPL` from this
        object.

        This is a bare-bones mechanism of creating a plot window from this
        object, which can then be moved around, zoomed, and on and on.  All
        behavior of the plot window is relegated to that routine.
        """
        pw = self._get_pw(fields, center, width, origin, axes_unit, 'Slice')
        return pw

class YTCuttingPlaneBase(YTSelectionContainer2D):
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
    _plane = None
    _top_node = "/CuttingPlanes"
    _key_fields = YTSelectionContainer2D._key_fields + ['pz','pdz']
    _type_name = "cutting"
    _con_args = ('normal', 'center')
    _container_fields = ("px", "py", "pz", "pdx", "pdy", "pdz")

    def __init__(self, normal, center, pf = None,
                 north_vector = None, field_parameters = None):
        YTSelectionContainer2D.__init__(self, 4, pf, field_parameters)
        self._set_center(center)
        self.set_field_parameter('center',center)
        # Let's set up our plane equation
        # ax + by + cz + d = 0
        self.orienter = Orientation(normal, north_vector = north_vector)
        self._norm_vec = self.orienter.normal_vector
        self._d = -1.0 * np.dot(self._norm_vec, self.center)
        self._x_vec = self.orienter.unit_vectors[0]
        self._y_vec = self.orienter.unit_vectors[1]
        # First we try all three, see which has the best result:
        vecs = np.identity(3)
        self._rot_mat = np.array([self._x_vec,self._y_vec,self._norm_vec])
        self._inv_mat = np.linalg.pinv(self._rot_mat)
        self.set_field_parameter('cp_x_vec',self._x_vec)
        self.set_field_parameter('cp_y_vec',self._y_vec)
        self.set_field_parameter('cp_z_vec',self._norm_vec)

    @property
    def normal(self):
        return self._norm_vec

    def to_frb(self, width, resolution, height=None):
        r"""This function returns an ObliqueFixedResolutionBuffer generated
        from this object.

        An ObliqueFixedResolutionBuffer is an object that accepts a
        variable-resolution 2D object and transforms it into an NxM bitmap that
        can be plotted, examined or processed.  This is a convenience function
        to return an FRB directly from an existing 2D data object.  Unlike the
        corresponding to_frb function for other AMR2DData objects, this does
        not accept a 'center' parameter as it is assumed to be centered at the
        center of the cutting plane.

        Parameters
        ----------
        width : width specifier
            This can either be a floating point value, in the native domain
            units of the simulation, or a tuple of the (value, unit) style.
            This will be the width of the FRB.
        height : height specifier, optional
            This will be the height of the FRB, by default it is equal to width.
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
        >>> write_image(np.log10(frb["Density"]), 'density_1pc.png')
        """
        if iterable(width):
            w, u = width
            width = w/self.pf[u]
        if height is None:
            height = width
        elif iterable(height):
            h, u = height
            height = h/self.pf[u]
        if not iterable(resolution):
            resolution = (resolution, resolution)
        from yt.visualization.fixed_resolution import ObliqueFixedResolutionBuffer
        bounds = (-width/2.0, width/2.0, -height/2.0, height/2.0)
        frb = ObliqueFixedResolutionBuffer(self, bounds, resolution)
        return frb

    def _generate_container_field(self, field):
        if self._current_chunk is None:
            self.hierarchy._identify_base_chunk(self)
        if field == "px":
            x = self._current_chunk.fcoords[:,0] - self.center[0]
            y = self._current_chunk.fcoords[:,1] - self.center[1]
            z = self._current_chunk.fcoords[:,2] - self.center[2]
            tr = np.zeros(self.size, dtype='float64')
            tr += x * self._x_vec[0]
            tr += y * self._x_vec[1]
            tr += z * self._x_vec[2]
            return tr
        elif field == "py":
            x = self._current_chunk.fcoords[:,0] - self.center[0]
            y = self._current_chunk.fcoords[:,1] - self.center[1]
            z = self._current_chunk.fcoords[:,2] - self.center[2]
            tr = np.zeros(self.size, dtype='float64')
            tr += x * self._y_vec[0]
            tr += y * self._y_vec[1]
            tr += z * self._y_vec[2]
            return tr
        elif field == "pz":
            x = self._current_chunk.fcoords[:,0] - self.center[0]
            y = self._current_chunk.fcoords[:,1] - self.center[1]
            z = self._current_chunk.fcoords[:,2] - self.center[2]
            tr = np.zeros(self.size, dtype='float64')
            tr += x * self._norm_vec[0]
            tr += y * self._norm_vec[1]
            tr += z * self._norm_vec[2]
            return tr
        elif field == "pdx":
            return self._current_chunk.fwidth[:,0] * 0.5
        elif field == "pdy":
            return self._current_chunk.fwidth[:,1] * 0.5
        elif field == "pdz":
            return self._current_chunk.fwidth[:,2] * 0.5
        else:
            raise KeyError(field)

    def to_pw(self, fields=None, center='c', width=None, axes_unit=None):
        r"""Create a :class:`~yt.visualization.plot_window.PWViewerMPL` from this
        object.

        This is a bare-bones mechanism of creating a plot window from this
        object, which can then be moved around, zoomed, and on and on.  All
        behavior of the plot window is relegated to that routine.
        """
        normal = self.normal
        center = self.center
        self.fields = [k for k in self.field_data.keys()
                       if k not in self._key_fields]
        from yt.visualization.plot_window import \
            GetObliqueWindowParameters, PWViewerMPL
        from yt.visualization.fixed_resolution import \
            ObliqueFixedResolutionBuffer
        (bounds, center_rot, units) = \
          GetObliqueWindowParameters(normal, center, width, self.pf)
        if axes_unit is None and units != ('1', '1'):
            axes_units = units
        pw = PWViewerMPL(
            self, bounds, fields=self.fields, origin='center-window',
            periodic=False, oblique=True,
            frb_generator=ObliqueFixedResolutionBuffer,
            plot_type='OffAxisSlice')
        pw.set_axes_unit(axes_unit)
        return pw

    def to_frb(self, width, resolution, height=None):
        r"""This function returns an ObliqueFixedResolutionBuffer generated
        from this object.

        An ObliqueFixedResolutionBuffer is an object that accepts a
        variable-resolution 2D object and transforms it into an NxM bitmap that
        can be plotted, examined or processed.  This is a convenience function
        to return an FRB directly from an existing 2D data object.  Unlike the
        corresponding to_frb function for other AMR2DData objects, this does
        not accept a 'center' parameter as it is assumed to be centered at the
        center of the cutting plane.

        Parameters
        ----------
        width : width specifier
            This can either be a floating point value, in the native domain
            units of the simulation, or a tuple of the (value, unit) style.
            This will be the width of the FRB.
        height : height specifier, optional
            This will be the height of the FRB, by default it is equal to width.
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
        >>> write_image(np.log10(frb["Density"]), 'density_1pc.png')
        """
        if iterable(width):
            w, u = width
            width = w/self.pf[u]
        if height is None:
            height = width
        elif iterable(height):
            h, u = height
            height = h/self.pf[u]
        if not iterable(resolution):
            resolution = (resolution, resolution)
        from yt.visualization.fixed_resolution import ObliqueFixedResolutionBuffer
        bounds = (-width/2.0, width/2.0, -height/2.0, height/2.0)
        frb = ObliqueFixedResolutionBuffer(self, bounds, resolution)
        return frb

class YTDiskBase(YTSelectionContainer3D):
    """
    By providing a *center*, a *normal*, a *radius* and a *height* we
    can define a cylinder of any proportion.  Only cells whose centers are
    within the cylinder will be selected.
    """
    _type_name = "disk"
    _con_args = ('center', '_norm_vec', '_radius', '_height')
    def __init__(self, center, normal, radius, height, fields=None,
                 pf=None, **kwargs):
        YTSelectionContainer3D.__init__(self, center, fields, pf, **kwargs)
        self._norm_vec = np.array(normal)/np.sqrt(np.dot(normal,normal))
        self.set_field_parameter("normal", self._norm_vec)
        self._height = fix_length(height, self.pf)
        self._radius = fix_length(radius, self.pf)
        self._d = -1.0 * np.dot(self._norm_vec, self.center)


class YTRegionBase(YTSelectionContainer3D):
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
    _type_name = "region"
    _con_args = ('center', 'left_edge', 'right_edge')
    def __init__(self, center, left_edge, right_edge, fields = None,
                 pf = None, **kwargs):
        YTSelectionContainer3D.__init__(self, center, fields, pf, **kwargs)
        self.left_edge = left_edge
        self.right_edge = right_edge

class YTDataCollectionBase(YTSelectionContainer3D):
    """
    By selecting an arbitrary *object_list*, we can act on those grids.
    Child cells are not returned.
    """
    _type_name = "data_collection"
    _con_args = ("_obj_list",)
    def __init__(self, center, obj_list, pf = None, field_parameters = None):
        YTSelectionContainer3D.__init__(self, center, pf, field_parameters)
        self._obj_ids = np.array([o.id - o._id_offset for o in obj_list],
                                dtype="int64")
        self._obj_list = obj_list

class YTSphereBase(YTSelectionContainer3D):
    """
    A sphere f points defined by a *center* and a *radius*.

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
    _type_name = "sphere"
    _con_args = ('center', 'radius')
    def __init__(self, center, radius, pf = None, field_parameters = None):
        super(YTSphereBase, self).__init__(center, pf, field_parameters)
        # Unpack the radius, if necessary
        radius = fix_length(radius, self.pf)
        #if radius < self.hierarchy.get_smallest_dx():
        #    raise YTSphereTooSmall(pf, radius, self.hierarchy.get_smallest_dx())
        self.set_field_parameter('radius',radius)
        self.radius = radius

class YTEllipsoidBase(YTSelectionContainer3D):
    """
    By providing a *center*,*A*,*B*,*C*,*e0*,*tilt* we
    can define a ellipsoid of any proportion.  Only cells whose
    centers are within the ellipsoid will be selected.

    Parameters
    ----------
    center : array_like
        The center of the ellipsoid.
    A : float
        The magnitude of the largest semi-major axis of the ellipsoid.
    B : float
        The magnitude of the medium semi-major axis of the ellipsoid.
    C : float
        The magnitude of the smallest semi-major axis of the ellipsoid.
    e0 : array_like (automatically normalized)
        the direction of the largest semi-major axis of the ellipsoid
    tilt : float
        After the rotation about the z-axis to allign e0 to x in the x-y
        plane, and then rotating about the y-axis to align e0 completely
        to the x-axis, tilt is the angle in radians remaining to
        rotate about the x-axis to align both e1 to the y-axis and e2 to
        the z-axis.
    Examples
    --------
    >>> pf = load("DD####/DD####")
    >>> c = [0.5,0.5,0.5]
    >>> ell = pf.h.ellipsoid(c, 0.1, 0.1, 0.1, np.array([0.1, 0.1, 0.1]), 0.2)
    """
    _type_name = "ellipsoid"
    _con_args = ('center', '_A', '_B', '_C', '_e0', '_tilt')
    def __init__(self, center, A, B, C, e0, tilt, fields=None,
                 pf=None, field_parameters = None):
        YTSelectionContainer3D.__init__(self, np.array(center), pf,
                                        field_parameters)
        # make sure the magnitudes of semi-major axes are in order
        if A<B or B<C:
            raise YTEllipsoidOrdering(pf, A, B, C)
        # make sure the smallest side is not smaller than dx
        #if C < self.hierarchy.get_smallest_dx():
        #    raise YTSphereTooSmall(pf, C, self.hierarchy.get_smallest_dx())
        self._A = A
        self._B = B
        self._C = C
        self._e0 = e0 = e0 / (e0**2.0).sum()**0.5
        self._tilt = tilt
        
        # find the t1 angle needed to rotate about z axis to align e0 to x
        t1 = np.arctan(e0[1] / e0[0])
        # rotate e0 by -t1
        RZ = get_rotation_matrix(t1, (0,0,1)).transpose()
        r1 = (e0 * RZ).sum(axis = 1)
        # find the t2 angle needed to rotate about y axis to align e0 to x
        t2 = np.arctan(-r1[2] / r1[0])
        """
        calculate the original e1
        given the tilt about the x axis when e0 was aligned 
        to x after t1, t2 rotations about z, y
        """
        RX = get_rotation_matrix(-tilt, (1, 0, 0)).transpose()
        RY = get_rotation_matrix(-t2,   (0, 1, 0)).transpose()
        RZ = get_rotation_matrix(-t1,   (0, 0, 1)).transpose()
        e1 = ((0, 1, 0) * RX).sum(axis=1)
        e1 = (e1 * RY).sum(axis=1)
        e1 = (e1 * RZ).sum(axis=1)
        e2 = np.cross(e0, e1)

        self._e1 = e1
        self._e2 = e2

        self.set_field_parameter('A', A)
        self.set_field_parameter('B', B)
        self.set_field_parameter('C', C)
        self.set_field_parameter('e0', e0)
        self.set_field_parameter('e1', e1)
        self.set_field_parameter('e2', e2)
