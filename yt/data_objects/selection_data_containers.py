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

import numpy as np

from yt.data_objects.data_containers import \
    YTSelectionContainer0D, YTSelectionContainer1D, \
    YTSelectionContainer2D, YTSelectionContainer3D
from yt.funcs import \
    ensure_list, \
    iterable, \
    validate_width_tuple, \
    fix_length
from yt.geometry.selection_routines import \
    points_in_cells
from yt.units.yt_array import \
    YTArray
from yt.utilities.exceptions import \
    YTSphereTooSmall, \
    YTIllDefinedCutRegion, \
    YTEllipsoidOrdering
from yt.utilities.minimal_representation import \
    MinimalSliceData
from yt.utilities.math_utils import get_rotation_matrix
from yt.utilities.orientation import Orientation


class YTPoint(YTSelectionContainer0D):
    """
    A 0-dimensional object defined by a single point

    Parameters
    ----------
    p: array_like
        A points defined within the domain.  If the domain is
        periodic its position will be corrected to lie inside
        the range [DLE,DRE) to ensure one and only one cell may
        match that point
    ds: Dataset, optional
        An optional dataset to use rather than self.ds
    field_parameters : dictionary
        A dictionary of field parameters than can be accessed by derived
        fields.
    data_source: optional
        Draw the selection from the provided data source rather than
        all data associated with the data_set

    Examples
    --------

    >>> import yt
    >>> ds = yt.load("RedshiftOutput0005")
    >>> c = [0.5,0.5,0.5]
    >>> point = ds.point(c)
    """
    _type_name = "point"
    _con_args = ('p',)
    def __init__(self, p, ds=None, field_parameters=None, data_source=None):
        super(YTPoint, self).__init__(ds, field_parameters, data_source)
        if isinstance(p, YTArray):
            # we pass p through ds.arr to ensure code units are attached
            self.p = self.ds.arr(p)
        else:
            self.p = self.ds.arr(p, 'code_length')

class YTOrthoRay(YTSelectionContainer1D):
    """
    This is an orthogonal ray cast through the entire domain, at a specific
    coordinate.

    This object is typically accessed through the `ortho_ray` object that
    hangs off of index objects.  The resulting arrays have their
    dimensionality reduced to one, and an ordered list of points at an
    (x,y) tuple along `axis` are available.

    Parameters
    ----------
    axis : int
        The axis along which to cast the ray.  Can be 0, 1, or 2 for x, y, z.
    coords : tuple of floats
        The (plane_x, plane_y) coordinates at which to cast the ray.  Note
        that this is in the plane coordinates: so if you are casting along
        x, this will be (y, z).  If you are casting along y, this will be
        (z, x).  If you are casting along z, this will be (x, y).
    ds: Dataset, optional
        An optional dataset to use rather than self.ds
    field_parameters : dictionary
         A dictionary of field parameters than can be accessed by derived
         fields.
    data_source: optional
        Draw the selection from the provided data source rather than
        all data associated with the data_set

    Examples
    --------

    >>> import yt
    >>> ds = yt.load("RedshiftOutput0005")
    >>> oray = ds.ortho_ray(0, (0.2, 0.74))
    >>> print oray["Density"]

    Note: The low-level data representation for rays are not guaranteed to be 
    spatially ordered.  In particular, with AMR datasets, higher resolution 
    data is tagged on to the end of the ray.  If you want this data 
    represented in a spatially ordered manner, manually sort it by the "t" 
    field, which is the value of the parametric variable that goes from 0 at 
    the start of the ray to 1 at the end:

    >>> my_ray = ds.ortho_ray(...)
    >>> ray_sort = np.argsort(my_ray["t"])
    >>> density = my_ray["density"][ray_sort]
    """
    _key_fields = ['x','y','z','dx','dy','dz']
    _type_name = "ortho_ray"
    _con_args = ('axis', 'coords')
    def __init__(self, axis, coords, ds=None, 
                 field_parameters=None, data_source=None):
        super(YTOrthoRay, self).__init__(ds, field_parameters, data_source)
        self.axis = axis
        xax = self.ds.coordinates.x_axis[self.axis]
        yax = self.ds.coordinates.y_axis[self.axis]
        self.px_ax = xax
        self.py_ax = yax
        # Even though we may not be using x,y,z we use them here.
        self.px_dx = 'd%s'%('xyz'[self.px_ax])
        self.py_dx = 'd%s'%('xyz'[self.py_ax])
        self.px, self.py = coords
        self.sort_by = 'xyz'[self.axis]

    @property
    def coords(self):
        return (self.px, self.py)

class YTRay(YTSelectionContainer1D):
    """
    This is an arbitrarily-aligned ray cast through the entire domain, at a
    specific coordinate.

    This object is typically accessed through the `ray` object that hangs
    off of index objects.  The resulting arrays have their
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
    ds: Dataset, optional
        An optional dataset to use rather than self.ds
    field_parameters : dictionary
         A dictionary of field parameters than can be accessed by derived
         fields.
    data_source: optional
        Draw the selection from the provided data source rather than
        all data associated with the data_set

    Examples
    --------

    >>> import yt
    >>> ds = yt.load("RedshiftOutput0005")
    >>> ray = ds.ray((0.2, 0.74, 0.11), (0.4, 0.91, 0.31))
    >>> print ray["Density"], ray["t"], ray["dts"]

    Note: The low-level data representation for rays are not guaranteed to be 
    spatially ordered.  In particular, with AMR datasets, higher resolution 
    data is tagged on to the end of the ray.  If you want this data 
    represented in a spatially ordered manner, manually sort it by the "t" 
    field, which is the value of the parametric variable that goes from 0 at 
    the start of the ray to 1 at the end:

    >>> my_ray = ds.ray(...)
    >>> ray_sort = np.argsort(my_ray["t"])
    >>> density = my_ray["density"][ray_sort]

"""
    _type_name = "ray"
    _con_args = ('start_point', 'end_point')
    _container_fields = ("t", "dts")
    def __init__(self, start_point, end_point, ds=None,
                 field_parameters=None, data_source=None):
        super(YTRay, self).__init__(ds, field_parameters, data_source)
        if isinstance(start_point, YTArray):
            self.start_point = \
              self.ds.arr(start_point).to("code_length")
        else:
            self.start_point = \
              self.ds.arr(start_point, 'code_length',
                          dtype='float64')
        if isinstance(end_point, YTArray):
            self.end_point = \
              self.ds.arr(end_point).to("code_length")
        else:
            self.end_point = \
              self.ds.arr(end_point, 'code_length',
                          dtype='float64')
        self.vec = self.end_point - self.start_point
        self._set_center(self.start_point)
        self.set_field_parameter('center', self.start_point)
        self._dts, self._ts = None, None

    def _generate_container_field(self, field):
        if self._current_chunk is None:
            self.index._identify_base_chunk(self)
        if field == "dts":
            return self._current_chunk.dtcoords
        elif field == "t":
            return self._current_chunk.tcoords
        else:
            raise KeyError(field)

class YTSlice(YTSelectionContainer2D):
    """
    This is a data object corresponding to a slice through the simulation
    domain.

    This object is typically accessed through the `slice` object that hangs
    off of index objects.  Slice is an orthogonal slice through the
    data, taking all the points at the finest resolution available and then
    indexing them.  It is more appropriately thought of as a slice
    'operator' than an object, however, as its field and coordinate can
    both change.

    Parameters
    ----------
    axis : int or char
        The axis along which to slice.  Can be 0, 1, or 2 for x, y, z.
    coord : float
        The coordinate along the axis at which to slice.  This is in
        "domain" coordinates.
    center : array_like, optional
        The 'center' supplied to fields that use it.  Note that this does
        not have to have `coord` as one value.  optional.
    ds: Dataset, optional
        An optional dataset to use rather than self.ds
    field_parameters : dictionary
         A dictionary of field parameters than can be accessed by derived
         fields.
    data_source: optional
        Draw the selection from the provided data source rather than
        all data associated with the data_set

    Examples
    --------

    >>> import yt
    >>> ds = yt.load("RedshiftOutput0005")
    >>> slice = ds.slice(0, 0.25)
    >>> print slice["Density"]
    """
    _top_node = "/Slices"
    _type_name = "slice"
    _con_args = ('axis', 'coord')
    _container_fields = ("px", "py", "pdx", "pdy")
    def __init__(self, axis, coord, center=None, ds=None,
                 field_parameters=None, data_source=None):
        YTSelectionContainer2D.__init__(self, axis, ds,
                                        field_parameters, data_source)
        self._set_center(center)
        self.coord = coord

    def _generate_container_field(self, field):
        xax = self.ds.coordinates.x_axis[self.axis]
        yax = self.ds.coordinates.y_axis[self.axis]
        if self._current_chunk is None:
            self.index._identify_base_chunk(self)
        if field == "px":
            return self._current_chunk.fcoords[:,xax]
        elif field == "py":
            return self._current_chunk.fcoords[:,yax]
        elif field == "pdx":
            return self._current_chunk.fwidth[:,xax] * 0.5
        elif field == "pdy":
            return self._current_chunk.fwidth[:,yax] * 0.5
        else:
            raise KeyError(field)

    @property
    def _mrep(self):
        return MinimalSliceData(self)

    def hub_upload(self):
        self._mrep.upload()

    def to_pw(self, fields=None, center='c', width=None, origin='center-window'):
        r"""Create a :class:`~yt.visualization.plot_window.PWViewerMPL` from this
        object.

        This is a bare-bones mechanism of creating a plot window from this
        object, which can then be moved around, zoomed, and on and on.  All
        behavior of the plot window is relegated to that routine.
        """
        pw = self._get_pw(fields, center, width, origin, 'Slice')
        return pw

    def plot(self, fields=None):
        if hasattr(self._data_source, "left_edge") and \
            hasattr(self._data_source, "right_edge"):
            left_edge = self._data_source.left_edge
            right_edge = self._data_source.right_edge
            center = (left_edge + right_edge)/2.0
            width = right_edge - left_edge
            xax = self.ds.coordinates.x_axis[self.axis]
            yax = self.ds.coordinates.y_axis[self.axis]
            lx, rx = left_edge[xax], right_edge[xax]
            ly, ry = left_edge[yax], right_edge[yax]
            width = (rx-lx), (ry-ly)
        else:
            width = self.ds.domain_width
            center = self.ds.domain_center
        pw = self._get_pw(fields, center, width, 'native', 'Slice')
        pw.show()
        return pw

class YTCuttingPlane(YTSelectionContainer2D):
    """
    This is a data object corresponding to an oblique slice through the
    simulation domain.

    This object is typically accessed through the `cutting` object
    that hangs off of index objects.  A cutting plane is an oblique
    plane through the data, defined by a normal vector and a coordinate.
    It attempts to guess an 'north' vector, which can be overridden, and
    then it pixelizes the appropriate data onto the plane without
    interpolation.

    Parameters
    ----------
    normal : array_like
        The vector that defines the desired plane.  For instance, the
        angular momentum of a sphere.
    center : array_like
        The center of the cutting plane, where the normal vector is anchored.
    north_vector: array_like, optional
        An optional vector to describe the north-facing direction in the resulting
        plane.
    ds: Dataset, optional
        An optional dataset to use rather than self.ds
    field_parameters : dictionary
         A dictionary of field parameters than can be accessed by derived
         fields.
    data_source: optional
        Draw the selection from the provided data source rather than
        all data associated with the data_set

    Notes
    -----

    This data object in particular can be somewhat expensive to create.
    It's also important to note that unlike the other 2D data objects, this
    object provides px, py, pz, as some cells may have a height from the
    plane.

    Examples
    --------

    >>> import yt
    >>> ds = yt.load("RedshiftOutput0005")
    >>> cp = ds.cutting([0.1, 0.2, -0.9], [0.5, 0.42, 0.6])
    >>> print cp["Density"]
    """
    _plane = None
    _top_node = "/CuttingPlanes"
    _key_fields = YTSelectionContainer2D._key_fields + ['pz','pdz']
    _type_name = "cutting"
    _con_args = ('normal', 'center')
    _tds_attrs = ("_inv_mat",)
    _tds_fields = ("x", "y", "z", "dx")
    _container_fields = ("px", "py", "pz", "pdx", "pdy", "pdz")
    def __init__(self, normal, center, north_vector=None,
                 ds=None, field_parameters=None, data_source=None):
        YTSelectionContainer2D.__init__(self, 4, ds,
                                        field_parameters, data_source)
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
        self._rot_mat = np.array([self._x_vec,self._y_vec,self._norm_vec])
        self._inv_mat = np.linalg.pinv(self._rot_mat)
        self.set_field_parameter('cp_x_vec',self._x_vec)
        self.set_field_parameter('cp_y_vec',self._y_vec)
        self.set_field_parameter('cp_z_vec',self._norm_vec)

    @property
    def normal(self):
        return self._norm_vec

    def _generate_container_field(self, field):
        if self._current_chunk is None:
            self.index._identify_base_chunk(self)
        if field == "px":
            x = self._current_chunk.fcoords[:,0] - self.center[0]
            y = self._current_chunk.fcoords[:,1] - self.center[1]
            z = self._current_chunk.fcoords[:,2] - self.center[2]
            tr = np.zeros(x.size, dtype='float64')
            tr = self.ds.arr(tr, "code_length")
            tr += x * self._x_vec[0]
            tr += y * self._x_vec[1]
            tr += z * self._x_vec[2]
            return tr
        elif field == "py":
            x = self._current_chunk.fcoords[:,0] - self.center[0]
            y = self._current_chunk.fcoords[:,1] - self.center[1]
            z = self._current_chunk.fcoords[:,2] - self.center[2]
            tr = np.zeros(x.size, dtype='float64')
            tr = self.ds.arr(tr, "code_length")
            tr += x * self._y_vec[0]
            tr += y * self._y_vec[1]
            tr += z * self._y_vec[2]
            return tr
        elif field == "pz":
            x = self._current_chunk.fcoords[:,0] - self.center[0]
            y = self._current_chunk.fcoords[:,1] - self.center[1]
            z = self._current_chunk.fcoords[:,2] - self.center[2]
            tr = np.zeros(x.size, dtype='float64')
            tr = self.ds.arr(tr, "code_length")
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
        self.fields = ensure_list(fields) + [k for k in self.field_data.keys()
                                             if k not in self._key_fields]
        from yt.visualization.plot_window import get_oblique_window_parameters, PWViewerMPL
        from yt.visualization.fixed_resolution import ObliqueFixedResolutionBuffer
        (bounds, center_rot) = get_oblique_window_parameters(normal, center, width, self.ds)
        pw = PWViewerMPL(
            self, bounds, fields=self.fields, origin='center-window', 
            periodic=False, oblique=True,
            frb_generator=ObliqueFixedResolutionBuffer, 
            plot_type='OffAxisSlice')
        if axes_unit is not None:
            pw.set_axes_unit(axes_unit)
        pw._setup_plots()
        return pw

    def to_frb(self, width, resolution, height=None, periodic=False):
        r"""This function returns an ObliqueFixedResolutionBuffer generated
        from this object.

        An ObliqueFixedResolutionBuffer is an object that accepts a
        variable-resolution 2D object and transforms it into an NxM bitmap that
        can be plotted, examined or processed.  This is a convenience function
        to return an FRB directly from an existing 2D data object.  Unlike the
        corresponding to_frb function for other YTSelectionContainer2D objects, 
        this does not accept a 'center' parameter as it is assumed to be 
        centered at the center of the cutting plane.

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
        periodic : boolean
            This can be true or false, and governs whether the pixelization
            will span the domain boundaries.

        Returns
        -------
        frb : :class:`~yt.visualization.fixed_resolution.ObliqueFixedResolutionBuffer`
            A fixed resolution buffer, which can be queried for fields.

        Examples
        --------

        >>> v, c = ds.find_max("density")
        >>> sp = ds.sphere(c, (100.0, 'au'))
        >>> L = sp.quantities.angular_momentum_vector()
        >>> cutting = ds.cutting(L, c)
        >>> frb = cutting.to_frb( (1.0, 'pc'), 1024)
        >>> write_image(np.log10(frb["Density"]), 'density_1pc.png')
        """
        if iterable(width):
            validate_width_tuple(width)
            width = self.ds.quan(width[0], width[1])
        if height is None:
            height = width
        elif iterable(height):
            validate_width_tuple(height)
            height = self.ds.quan(height[0], height[1])
        if not iterable(resolution):
            resolution = (resolution, resolution)
        from yt.visualization.fixed_resolution import ObliqueFixedResolutionBuffer
        bounds = (-width/2.0, width/2.0, -height/2.0, height/2.0)
        frb = ObliqueFixedResolutionBuffer(self, bounds, resolution,
                                           periodic=periodic)
        return frb

class YTDisk(YTSelectionContainer3D):
    """
    By providing a *center*, a *normal*, a *radius* and a *height* we
    can define a cylinder of any proportion.  Only cells whose centers are
    within the cylinder will be selected.

    Parameters
    ----------
    center : array_like
        coordinate to which the normal, radius, and height all reference
    normal : array_like
        the normal vector defining the direction of lengthwise part of the 
        cylinder
    radius : float
        the radius of the cylinder
    height : float
        the distance from the midplane of the cylinder to the top and 
        bottom planes
    fields : array of fields, optional
        any fields to be pre-loaded in the cylinder object
    ds: Dataset, optional
        An optional dataset to use rather than self.ds
    field_parameters : dictionary
         A dictionary of field parameters than can be accessed by derived
         fields.
    data_source: optional
        Draw the selection from the provided data source rather than
        all data associated with the data_set

    Examples
    --------

    >>> import yt
    >>> ds = yt.load("RedshiftOutput0005")
    >>> c = [0.5,0.5,0.5]
    >>> disk = ds.disk(c, [1,0,0], (1, 'kpc'), (10, 'kpc'))
    """
    _type_name = "disk"
    _con_args = ('center', '_norm_vec', 'radius', 'height')
    def __init__(self, center, normal, radius, height, fields=None,
                 ds=None, field_parameters=None, data_source=None):
        YTSelectionContainer3D.__init__(self, center, ds,
                                        field_parameters, data_source)
        self._norm_vec = np.array(normal)/np.sqrt(np.dot(normal,normal))
        self.set_field_parameter("normal", self._norm_vec)
        self.set_field_parameter("center", self.center)
        self.height = fix_length(height, self.ds)
        self.radius = fix_length(radius, self.ds)
        self._d = -1.0 * np.dot(self._norm_vec, self.center)

class YTRegion(YTSelectionContainer3D):
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
    def __init__(self, center, left_edge, right_edge, fields=None,
                 ds=None, field_parameters=None, data_source=None):
        YTSelectionContainer3D.__init__(self, center, ds,
                                        field_parameters, data_source)
        if not isinstance(left_edge, YTArray):
            self.left_edge = self.ds.arr(left_edge, 'code_length')
        else:
            self.left_edge = left_edge
        if not isinstance(right_edge, YTArray):
            self.right_edge = self.ds.arr(right_edge, 'code_length')
        else:
            self.right_edge = right_edge

class YTDataCollection(YTSelectionContainer3D):
    """
    By selecting an arbitrary *object_list*, we can act on those grids.
    Child cells are not returned.
    """
    _type_name = "data_collection"
    _con_args = ("_obj_list",)
    def __init__(self, obj_list, ds=None, field_parameters=None,
                 data_source=None, center=None):
        YTSelectionContainer3D.__init__(self, center, ds,
                                        field_parameters, data_source)
        self._obj_ids = np.array([o.id - o._id_offset for o in obj_list],
                                dtype="int64")
        self._obj_list = obj_list

class YTSphere(YTSelectionContainer3D):
    """
    A sphere of points defined by a *center* and a *radius*.

    Parameters
    ----------
    center : array_like
        The center of the sphere.
    radius : float
        The radius of the sphere.

    Examples
    --------

    >>> import yt
    >>> ds = yt.load("RedshiftOutput0005")
    >>> c = [0.5,0.5,0.5]
    >>> sphere = ds.sphere(c, (1., "kpc"))
    """
    _type_name = "sphere"
    _con_args = ('center', 'radius')
    def __init__(self, center, radius, ds=None,
                 field_parameters=None, data_source=None):
        super(YTSphere, self).__init__(center, ds,
                                           field_parameters, data_source)
        # Unpack the radius, if necessary
        radius = fix_length(radius, self.ds)
        if radius < self.index.get_smallest_dx():
            raise YTSphereTooSmall(ds, radius.in_units("code_length"),
                                   self.index.get_smallest_dx().in_units("code_length"))
        self.set_field_parameter('radius',radius)
        self.set_field_parameter("center", self.center)
        self.radius = radius

class YTEllipsoid(YTSelectionContainer3D):
    """
    By providing a *center*,*A*,*B*,*C*,*e0*,*tilt* we
    can define a ellipsoid of any proportion.  Only cells whose
    centers are within the ellipsoid will be selected.

    Parameters
    ----------
    center : array_like
        The center of the ellipsoid.
    A : float
        The magnitude of the largest axis (semi-major) of the ellipsoid.
    B : float
        The magnitude of the medium axis (semi-medium) of the ellipsoid.
    C : float
        The magnitude of the smallest axis (semi-minor) of the ellipsoid.
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

    >>> import yt
    >>> ds = yt.load("RedshiftOutput0005")
    >>> c = [0.5,0.5,0.5]
    >>> ell = ds.ellipsoid(c, 0.1, 0.1, 0.1, np.array([0.1, 0.1, 0.1]), 0.2)
    """
    _type_name = "ellipsoid"
    _con_args = ('center', '_A', '_B', '_C', '_e0', '_tilt')
    def __init__(self, center, A, B, C, e0, tilt, fields=None,
                 ds=None, field_parameters=None, data_source=None):
        YTSelectionContainer3D.__init__(self, center, ds,
                                        field_parameters, data_source)
        # make sure the magnitudes of semi-major axes are in order
        if A<B or B<C:
            raise YTEllipsoidOrdering(ds, A, B, C)
        # make sure the smallest side is not smaller than dx
        self._A = self.ds.quan(A, 'code_length')
        self._B = self.ds.quan(B, 'code_length')
        self._C = self.ds.quan(C, 'code_length')
        if self._C < self.index.get_smallest_dx():
            raise YTSphereTooSmall(self.ds, self._C, self.index.get_smallest_dx())
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

class YTCutRegion(YTSelectionContainer3D):
    """
    This is a data object designed to allow individuals to apply logical
    operations to fields and filter as a result of those cuts.

    Parameters
    ----------
    data_source : YTSelectionContainer3D
        The object to which cuts will be applied.
    conditionals : list of strings
        A list of conditionals that will be evaluated.  In the namespace
        available, these conditionals will have access to 'obj' which is a data
        object of unknown shape, and they must generate a boolean array.  For
        instance, conditionals = ["obj['temperature'] < 1e3"]

    Examples
    --------

    >>> import yt
    >>> ds = yt.load("RedshiftOutput0005")
    >>> sp = ds.sphere("max", (1.0, 'Mpc'))
    >>> cr = ds.cut_region(sp, ["obj['temperature'] < 1e3"])
    """
    _type_name = "cut_region"
    _con_args = ("base_object", "conditionals")
    def __init__(self, data_source, conditionals, ds=None,
                 field_parameters=None, base_object=None):
        if base_object is not None:
            # passing base_object explicitly has been deprecated,
            # but we handle it here for backward compatibility
            if data_source is not None:
                raise RuntimeError(
                    "Cannot use both base_object and data_source")
            data_source=base_object
        super(YTCutRegion, self).__init__(
            data_source.center, ds, field_parameters, data_source=data_source)
        self.conditionals = ensure_list(conditionals)
        self.base_object = data_source
        self._selector = None
        # Need to interpose for __getitem__, fwidth, fcoords, icoords, iwidth,
        # ires and get_data

    def chunks(self, fields, chunking_style, **kwargs):
        # We actually want to chunk the sub-chunk, not ourselves.  We have no
        # chunks to speak of, as we do not data IO.
        for chunk in self.index._chunk(self.base_object,
                                       chunking_style,
                                       **kwargs):
            with self.base_object._chunked_read(chunk):
                with self._chunked_read(chunk):
                    self.get_data(fields)
                    yield self

    def get_data(self, fields = None):
        fields = ensure_list(fields)
        self.base_object.get_data(fields)
        ind = self._cond_ind
        for field in fields:
            f = self.base_object[field]
            if f.shape != ind.shape:
                parent = getattr(self, "parent", self.base_object)
                self.field_data[field] = parent[field][self._part_ind]
            else:
                self.field_data[field] = self.base_object[field][ind]

    @property
    def blocks(self):
        # We have to take a slightly different approach here.  Note that all
        # that .blocks has to yield is a 3D array and a mask.
        for obj, m in self.base_object.blocks:
            m = m.copy()
            with obj._field_parameter_state(self.field_parameters):
                for cond in self.conditionals:
                    ss = eval(cond)
                    m = np.logical_and(m, ss, m)
            if not np.any(m): continue
            yield obj, m

    @property
    def _cond_ind(self):
        ind = None
        obj = self.base_object
        with obj._field_parameter_state(self.field_parameters):
            for cond in self.conditionals:
                res = eval(cond)
                if ind is None: ind = res
                if ind.shape != res.shape:
                    raise YTIllDefinedCutRegion(self.conditionals)
                np.logical_and(res, ind, ind)
        return ind

    _particle_mask = None
    @property
    def _part_ind(self):
        if self._particle_mask is None:
            parent = getattr(self, "parent", self.base_object)
            units = "code_length"
            mask = points_in_cells(
                self["x"].to(units), self["y"].to(units),
                self["z"].to(units), self["dx"].to(units),
                self["dy"].to(units), self["dz"].to(units),
                parent["particle_position_x"].to(units),
                parent["particle_position_y"].to(units),
                parent["particle_position_z"].to(units))
            self._particle_mask = mask
        return self._particle_mask

    @property
    def icoords(self):
        return self.base_object.icoords[self._cond_ind,:]

    @property
    def fcoords(self):
        return self.base_object.fcoords[self._cond_ind,:]

    @property
    def ires(self):
        return self.base_object.ires[self._cond_ind]

    @property
    def fwidth(self):
        return self.base_object.fwidth[self._cond_ind,:]
