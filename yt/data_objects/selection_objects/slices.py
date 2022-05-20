import numpy as np

from yt.data_objects.selection_objects.data_selection_objects import (
    YTSelectionContainer,
    YTSelectionContainer2D,
)
from yt.data_objects.static_output import Dataset
from yt.funcs import (
    is_sequence,
    iter_fields,
    validate_3d_array,
    validate_axis,
    validate_center,
    validate_float,
    validate_object,
    validate_width_tuple,
)
from yt.utilities.exceptions import YTNotInsideNotebook
from yt.utilities.minimal_representation import MinimalSliceData
from yt.utilities.orientation import Orientation


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
    ds: ~yt.data_objects.static_output.Dataset, optional
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
    >>> print(slice[("gas", "density")])
    """

    _top_node = "/Slices"
    _type_name = "slice"
    _con_args = ("axis", "coord")
    _container_fields = ("px", "py", "pz", "pdx", "pdy", "pdz")

    def __init__(
        self, axis, coord, center=None, ds=None, field_parameters=None, data_source=None
    ):
        validate_axis(ds, axis)
        validate_float(coord)
        # center is an optional parameter
        if center is not None:
            validate_center(center)
        validate_object(ds, Dataset)
        validate_object(field_parameters, dict)
        validate_object(data_source, YTSelectionContainer)
        YTSelectionContainer2D.__init__(self, axis, ds, field_parameters, data_source)
        self._set_center(center)
        self.coord = coord

    def _generate_container_field(self, field):
        xax = self.ds.coordinates.x_axis[self.axis]
        yax = self.ds.coordinates.y_axis[self.axis]
        if self._current_chunk is None:
            self.index._identify_base_chunk(self)
        if field == "px":
            return self._current_chunk.fcoords[:, xax]
        elif field == "py":
            return self._current_chunk.fcoords[:, yax]
        elif field == "pz":
            return self._current_chunk.fcoords[:, self.axis]
        elif field == "pdx":
            return self._current_chunk.fwidth[:, xax] * 0.5
        elif field == "pdy":
            return self._current_chunk.fwidth[:, yax] * 0.5
        elif field == "pdz":
            return self._current_chunk.fwidth[:, self.axis] * 0.5
        else:
            raise KeyError(field)

    @property
    def _mrep(self):
        return MinimalSliceData(self)

    def to_pw(self, fields=None, center="c", width=None, origin="center-window"):
        r"""Create a :class:`~yt.visualization.plot_window.PWViewerMPL` from this
        object.

        This is a bare-bones mechanism of creating a plot window from this
        object, which can then be moved around, zoomed, and on and on.  All
        behavior of the plot window is relegated to that routine.
        """
        pw = self._get_pw(fields, center, width, origin, "Slice")
        return pw

    def plot(self, fields=None):
        if hasattr(self._data_source, "left_edge") and hasattr(
            self._data_source, "right_edge"
        ):
            left_edge = self._data_source.left_edge
            right_edge = self._data_source.right_edge
            center = (left_edge + right_edge) / 2.0
            width = right_edge - left_edge
            xax = self.ds.coordinates.x_axis[self.axis]
            yax = self.ds.coordinates.y_axis[self.axis]
            lx, rx = left_edge[xax], right_edge[xax]
            ly, ry = left_edge[yax], right_edge[yax]
            width = (rx - lx), (ry - ly)
        else:
            width = self.ds.domain_width
            center = self.ds.domain_center
        pw = self._get_pw(fields, center, width, "native", "Slice")
        try:
            pw.show()
        except YTNotInsideNotebook:
            pass
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
    ds: ~yt.data_objects.static_output.Dataset, optional
        An optional dataset to use rather than self.ds
    field_parameters : dictionary
         A dictionary of field parameters than can be accessed by derived
         fields.
    data_source: optional
        Draw the selection from the provided data source rather than
        all data associated with the dataset

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
    >>> print(cp[("gas", "density")])
    """

    _plane = None
    _top_node = "/CuttingPlanes"
    _key_fields = YTSelectionContainer2D._key_fields + ["pz", "pdz"]
    _type_name = "cutting"
    _con_args = ("normal", "center")
    _tds_attrs = ("_inv_mat",)
    _tds_fields = ("x", "y", "z", "dx")
    _container_fields = ("px", "py", "pz", "pdx", "pdy", "pdz")

    def __init__(
        self,
        normal,
        center,
        north_vector=None,
        ds=None,
        field_parameters=None,
        data_source=None,
    ):
        validate_3d_array(normal)
        validate_center(center)
        if north_vector is not None:
            validate_3d_array(north_vector)
        validate_object(ds, Dataset)
        validate_object(field_parameters, dict)
        validate_object(data_source, YTSelectionContainer)
        YTSelectionContainer2D.__init__(self, 4, ds, field_parameters, data_source)
        self._set_center(center)
        self.set_field_parameter("center", center)
        # Let's set up our plane equation
        # ax + by + cz + d = 0
        self.orienter = Orientation(normal, north_vector=north_vector)
        self._norm_vec = self.orienter.normal_vector
        self._d = -1.0 * np.dot(self._norm_vec, self.center)
        self._x_vec = self.orienter.unit_vectors[0]
        self._y_vec = self.orienter.unit_vectors[1]
        # First we try all three, see which has the best result:
        self._rot_mat = np.array([self._x_vec, self._y_vec, self._norm_vec])
        self._inv_mat = np.linalg.pinv(self._rot_mat)
        self.set_field_parameter("cp_x_vec", self._x_vec)
        self.set_field_parameter("cp_y_vec", self._y_vec)
        self.set_field_parameter("cp_z_vec", self._norm_vec)

    @property
    def normal(self):
        return self._norm_vec

    def _generate_container_field(self, field):
        if self._current_chunk is None:
            self.index._identify_base_chunk(self)
        if field == "px":
            x = self._current_chunk.fcoords[:, 0] - self.center[0]
            y = self._current_chunk.fcoords[:, 1] - self.center[1]
            z = self._current_chunk.fcoords[:, 2] - self.center[2]
            tr = np.zeros(x.size, dtype="float64")
            tr = self.ds.arr(tr, "code_length")
            tr += x * self._x_vec[0]
            tr += y * self._x_vec[1]
            tr += z * self._x_vec[2]
            return tr
        elif field == "py":
            x = self._current_chunk.fcoords[:, 0] - self.center[0]
            y = self._current_chunk.fcoords[:, 1] - self.center[1]
            z = self._current_chunk.fcoords[:, 2] - self.center[2]
            tr = np.zeros(x.size, dtype="float64")
            tr = self.ds.arr(tr, "code_length")
            tr += x * self._y_vec[0]
            tr += y * self._y_vec[1]
            tr += z * self._y_vec[2]
            return tr
        elif field == "pz":
            x = self._current_chunk.fcoords[:, 0] - self.center[0]
            y = self._current_chunk.fcoords[:, 1] - self.center[1]
            z = self._current_chunk.fcoords[:, 2] - self.center[2]
            tr = np.zeros(x.size, dtype="float64")
            tr = self.ds.arr(tr, "code_length")
            tr += x * self._norm_vec[0]
            tr += y * self._norm_vec[1]
            tr += z * self._norm_vec[2]
            return tr
        elif field == "pdx":
            return self._current_chunk.fwidth[:, 0] * 0.5
        elif field == "pdy":
            return self._current_chunk.fwidth[:, 1] * 0.5
        elif field == "pdz":
            return self._current_chunk.fwidth[:, 2] * 0.5
        else:
            raise KeyError(field)

    def to_pw(self, fields=None, center="c", width=None, axes_unit=None):
        r"""Create a :class:`~yt.visualization.plot_window.PWViewerMPL` from this
        object.

        This is a bare-bones mechanism of creating a plot window from this
        object, which can then be moved around, zoomed, and on and on.  All
        behavior of the plot window is relegated to that routine.
        """
        normal = self.normal
        center = self.center
        self.fields = list(iter_fields(fields)) + [
            k for k in self.field_data.keys() if k not in self._key_fields
        ]
        from yt.visualization.fixed_resolution import FixedResolutionBuffer
        from yt.visualization.plot_window import (
            PWViewerMPL,
            get_oblique_window_parameters,
        )

        (bounds, center_rot) = get_oblique_window_parameters(
            normal, center, width, self.ds
        )
        pw = PWViewerMPL(
            self,
            bounds,
            fields=self.fields,
            origin="center-window",
            periodic=False,
            oblique=True,
            frb_generator=FixedResolutionBuffer,
            plot_type="OffAxisSlice",
        )
        if axes_unit is not None:
            pw.set_axes_unit(axes_unit)
        pw._setup_plots()
        return pw

    def to_frb(self, width, resolution, height=None, periodic=False):
        r"""This function returns a FixedResolutionBuffer generated from this
        object.

        An FixedResolutionBuffer is an object that accepts a
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
        frb : :class:`~yt.visualization.fixed_resolution.FixedResolutionBuffer`
            A fixed resolution buffer, which can be queried for fields.

        Examples
        --------

        >>> v, c = ds.find_max(("gas", "density"))
        >>> sp = ds.sphere(c, (100.0, "au"))
        >>> L = sp.quantities.angular_momentum_vector()
        >>> cutting = ds.cutting(L, c)
        >>> frb = cutting.to_frb((1.0, "pc"), 1024)
        >>> write_image(np.log10(frb[("gas", "density")]), "density_1pc.png")
        """
        if is_sequence(width):
            validate_width_tuple(width)
            width = self.ds.quan(width[0], width[1])
        if height is None:
            height = width
        elif is_sequence(height):
            validate_width_tuple(height)
            height = self.ds.quan(height[0], height[1])
        if not is_sequence(resolution):
            resolution = (resolution, resolution)
        from yt.visualization.fixed_resolution import FixedResolutionBuffer

        bounds = (-width / 2.0, width / 2.0, -height / 2.0, height / 2.0)
        frb = FixedResolutionBuffer(self, bounds, resolution, periodic=periodic)
        return frb
