import abc
import weakref
from numbers import Number
from typing import Tuple

import numpy as np
from packaging.version import Version

from yt.funcs import fix_unitary, is_sequence, validate_width_tuple
from yt.units.yt_array import YTArray, YTQuantity
from yt.utilities.exceptions import YTCoordinateNotImplemented, YTInvalidWidthError
from yt.visualization._commons import MPL_VERSION


def _unknown_coord(field, data):
    raise YTCoordinateNotImplemented


def _get_coord_fields(axi, units="code_length"):
    def _dds(field, data):
        rv = data.ds.arr(data.fwidth[..., axi].copy(), units)
        return data._reshape_vals(rv)

    def _coords(field, data):
        rv = data.ds.arr(data.fcoords[..., axi].copy(), units)
        return data._reshape_vals(rv)

    return _dds, _coords


def _get_vert_fields(axi, units="code_length"):
    def _vert(field, data):
        rv = data.ds.arr(data.fcoords_vertex[..., axi].copy(), units)
        return rv

    return _vert


def _setup_dummy_cartesian_coords_and_widths(registry, axes: Tuple[str]):
    for ax in axes:
        registry.add_field(
            ("index", f"d{ax}"), sampling_type="cell", function=_unknown_coord
        )
        registry.add_field(("index", ax), sampling_type="cell", function=_unknown_coord)


def _setup_polar_coordinates(registry, axis_id):
    f1, f2 = _get_coord_fields(axis_id["r"])
    registry.add_field(
        ("index", "dr"),
        sampling_type="cell",
        function=f1,
        display_field=False,
        units="code_length",
    )

    registry.add_field(
        ("index", "r"),
        sampling_type="cell",
        function=f2,
        display_field=False,
        units="code_length",
    )

    f1, f2 = _get_coord_fields(axis_id["theta"], "dimensionless")
    registry.add_field(
        ("index", "dtheta"),
        sampling_type="cell",
        function=f1,
        display_field=False,
        units="dimensionless",
    )

    registry.add_field(
        ("index", "theta"),
        sampling_type="cell",
        function=f2,
        display_field=False,
        units="dimensionless",
    )

    def _path_r(field, data):
        return data["index", "dr"]

    registry.add_field(
        ("index", "path_element_r"),
        sampling_type="cell",
        function=_path_r,
        units="code_length",
    )

    def _path_theta(field, data):
        # Note: this already assumes cell-centered
        return data["index", "r"] * data["index", "dtheta"]

    registry.add_field(
        ("index", "path_element_theta"),
        sampling_type="cell",
        function=_path_theta,
        units="code_length",
    )


def validate_sequence_width(width, ds, unit=None):
    if isinstance(width[0], tuple) and isinstance(width[1], tuple):
        validate_width_tuple(width[0])
        validate_width_tuple(width[1])
        return (
            ds.quan(width[0][0], fix_unitary(width[0][1])),
            ds.quan(width[1][0], fix_unitary(width[1][1])),
        )
    elif isinstance(width[0], Number) and isinstance(width[1], Number):
        return (ds.quan(width[0], "code_length"), ds.quan(width[1], "code_length"))
    elif isinstance(width[0], YTQuantity) and isinstance(width[1], YTQuantity):
        return (ds.quan(width[0]), ds.quan(width[1]))
    else:
        validate_width_tuple(width)
        # If width and unit are both valid width tuples, we
        # assume width controls x and unit controls y
        try:
            validate_width_tuple(unit)
            return (
                ds.quan(width[0], fix_unitary(width[1])),
                ds.quan(unit[0], fix_unitary(unit[1])),
            )
        except YTInvalidWidthError:
            return (
                ds.quan(width[0], fix_unitary(width[1])),
                ds.quan(width[0], fix_unitary(width[1])),
            )


class CoordinateHandler(abc.ABC):
    name: str

    def __init__(self, ds, ordering):
        self.ds = weakref.proxy(ds)
        self.axis_order = ordering

    @abc.abstractmethod
    def setup_fields(self):
        # This should return field definitions for x, y, z, r, theta, phi
        pass

    @abc.abstractmethod
    def pixelize(self, dimension, data_source, field, bounds, size, antialias=True):
        # This should *actually* be a pixelize call, not just returning the
        # pixelizer
        pass

    @abc.abstractmethod
    def pixelize_line(self, field, start_point, end_point, npoints):
        pass

    def distance(self, start, end):
        p1 = self.convert_to_cartesian(start)
        p2 = self.convert_to_cartesian(end)
        return np.sqrt(((p1 - p2) ** 2.0).sum())

    @abc.abstractmethod
    def convert_from_cartesian(self, coord):
        pass

    @abc.abstractmethod
    def convert_to_cartesian(self, coord):
        pass

    @abc.abstractmethod
    def convert_to_cylindrical(self, coord):
        pass

    @abc.abstractmethod
    def convert_from_cylindrical(self, coord):
        pass

    @abc.abstractmethod
    def convert_to_spherical(self, coord):
        pass

    @abc.abstractmethod
    def convert_from_spherical(self, coord):
        pass

    _data_projection = None

    @property
    def data_projection(self):
        if self._data_projection is not None:
            return self._data_projection
        dpj = {}
        for ax in self.axis_order:
            dpj[ax] = None
        self._data_projection = dpj
        return dpj

    _data_transform = None

    @property
    def data_transform(self):
        if self._data_transform is not None:
            return self._data_transform
        dtx = {}
        for ax in self.axis_order:
            dtx[ax] = None
        self._data_transform = dtx
        return dtx

    _axis_name = None

    @property
    def axis_name(self):
        if self._axis_name is not None:
            return self._axis_name
        an = {}
        for axi, ax in enumerate(self.axis_order):
            an[axi] = ax
            an[ax] = ax
            an[ax.capitalize()] = ax
        self._axis_name = an
        return an

    _axis_id = None

    @property
    def axis_id(self):
        if self._axis_id is not None:
            return self._axis_id
        ai = {}
        for axi, ax in enumerate(self.axis_order):
            ai[ax] = ai[axi] = axi
        self._axis_id = ai
        return ai

    _image_axis_name = None

    @property
    def image_axis_name(self):
        # Default
        if self._image_axis_name is not None:
            return self._image_axis_name
        self._image_axis_name = rv = {}
        for i in range(3):
            rv[i] = (self.axis_name[self.x_axis[i]], self.axis_name[self.y_axis[i]])
            rv[self.axis_name[i]] = rv[i]
            rv[self.axis_name[i].capitalize()] = rv[i]
        return rv

    _x_axis = None

    @property
    def x_axis(self):
        if self._x_axis is not None:
            return self._x_axis
        ai = self.axis_id
        xa = {}
        for a1, a2 in self._x_pairs:
            xa[a1] = xa[ai[a1]] = ai[a2]
        self._x_axis = xa
        return xa

    _y_axis = None

    @property
    def y_axis(self):
        if self._y_axis is not None:
            return self._y_axis
        ai = self.axis_id
        ya = {}
        for a1, a2 in self._y_pairs:
            ya[a1] = ya[ai[a1]] = ai[a2]
        self._y_axis = ya
        return ya

    @property
    @abc.abstractproperty
    def period(self):
        pass

    def sanitize_depth(self, depth):
        if is_sequence(depth):
            validate_width_tuple(depth)
            depth = (self.ds.quan(depth[0], fix_unitary(depth[1])),)
        elif isinstance(depth, Number):
            depth = (
                self.ds.quan(depth, "code_length", registry=self.ds.unit_registry),
            )
        elif isinstance(depth, YTQuantity):
            depth = (depth,)
        else:
            raise YTInvalidWidthError(depth)
        return depth

    def sanitize_width(self, axis, width, depth):
        if width is None:
            # initialize the index if it is not already initialized
            self.ds.index
            # Default to code units
            if not is_sequence(axis):
                xax = self.x_axis[axis]
                yax = self.y_axis[axis]
                w = self.ds.domain_width[np.array([xax, yax])]
            else:
                # axis is actually the normal vector
                # for an off-axis data object.
                mi = np.argmin(self.ds.domain_width)
                w = self.ds.domain_width[np.array((mi, mi))]
            width = (w[0], w[1])
        elif is_sequence(width):
            width = validate_sequence_width(width, self.ds)
        elif isinstance(width, YTQuantity):
            width = (width, width)
        elif isinstance(width, Number):
            width = (
                self.ds.quan(width, "code_length"),
                self.ds.quan(width, "code_length"),
            )
        else:
            raise YTInvalidWidthError(width)
        if depth is not None:
            depth = self.sanitize_depth(depth)
            return width + depth
        return width

    def sanitize_center(self, center, axis):
        if isinstance(center, str):
            if center.lower() == "m" or center.lower() == "max":
                v, center = self.ds.find_max(("gas", "density"))
                center = self.ds.arr(center, "code_length")
            elif center.lower() == "c" or center.lower() == "center":
                # domain_left_edge and domain_right_edge might not be
                # initialized until we create the index, so create it
                self.ds.index
                center = (self.ds.domain_left_edge + self.ds.domain_right_edge) / 2
            else:
                raise RuntimeError(f'center keyword "{center}" not recognized')
        elif isinstance(center, YTArray):
            return self.ds.arr(center), self.convert_to_cartesian(center)
        elif is_sequence(center):
            if isinstance(center[0], str) and isinstance(center[1], str):
                if center[0].lower() == "min":
                    v, center = self.ds.find_min(center[1])
                elif center[0].lower() == "max":
                    v, center = self.ds.find_max(center[1])
                else:
                    raise RuntimeError(f'center keyword "{center}" not recognized')
                center = self.ds.arr(center, "code_length")
            elif is_sequence(center[0]) and isinstance(center[1], str):
                center = self.ds.arr(center[0], center[1])
            else:
                center = self.ds.arr(center, "code_length")
        else:
            raise RuntimeError(f'center keyword "{center}" not recognized')
        # This has to return both a center and a display_center
        display_center = self.convert_to_cartesian(center)
        return center, display_center

    def sanitize_buffer_fill_values(self, buff):
        """Replace nans with +inf in buff, if all valid values are positive"""
        # In buffer with only positive values, maplotlib will raise a warning
        # if nan is used as a filler, while it tolerates np.inf just fine
        # This hack is however not necessary since Matpltolib 3.2
        if MPL_VERSION >= Version("3.2"):
            return
        minval = buff[~np.isnan(buff)].min()
        if minval >= 0:
            buff[np.isnan(buff)] = np.inf


def cartesian_to_cylindrical(coord, center=(0, 0, 0)):
    c2 = np.zeros_like(coord)
    if not isinstance(center, YTArray):
        center = center * coord.uq
    c2[..., 0] = (
        (coord[..., 0] - center[0]) ** 2.0 + (coord[..., 1] - center[1]) ** 2.0
    ) ** 0.5
    c2[..., 1] = coord[..., 2]  # rzt
    c2[..., 2] = np.arctan2(coord[..., 1] - center[1], coord[..., 0] - center[0])
    return c2


def cylindrical_to_cartesian(coord, center=(0, 0, 0)):
    c2 = np.zeros_like(coord)
    if not isinstance(center, YTArray):
        center = center * coord.uq
    c2[..., 0] = np.cos(coord[..., 0]) * coord[..., 1] + center[0]
    c2[..., 1] = np.sin(coord[..., 0]) * coord[..., 1] + center[1]
    c2[..., 2] = coord[..., 2]
    return c2


def _get_polar_bounds(self: CoordinateHandler, axes: Tuple[str, str]):
    # a small helper function that is needed by two unrelated classes
    ri = self.axis_id[axes[0]]
    pi = self.axis_id[axes[1]]
    rmin = self.ds.domain_left_edge[ri]
    rmax = self.ds.domain_right_edge[ri]
    phimin = self.ds.domain_left_edge[pi]
    phimax = self.ds.domain_right_edge[pi]
    corners = [
        (rmin, phimin),
        (rmin, phimax),
        (rmax, phimin),
        (rmax, phimax),
    ]

    def to_polar_plane(r, phi):
        x = r * np.cos(phi)
        y = r * np.sin(phi)
        return x, y

    conic_corner_coords = [to_polar_plane(*corner) for corner in corners]

    phimin = phimin.d
    phimax = phimax.d

    if phimin <= np.pi <= phimax:
        xxmin = -rmax
    else:
        xxmin = min(xx for xx, yy in conic_corner_coords)

    if phimin <= 0 <= phimax:
        xxmax = rmax
    else:
        xxmax = max(xx for xx, yy in conic_corner_coords)

    if phimin <= 3 * np.pi / 2 <= phimax:
        yymin = -rmax
    else:
        yymin = min(yy for xx, yy in conic_corner_coords)

    if phimin <= np.pi / 2 <= phimax:
        yymax = rmax
    else:
        yymax = max(yy for xx, yy in conic_corner_coords)

    return xxmin, xxmax, yymin, yymax
