"""
Coordinate handler base class.




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
import weakref
from numbers import Number

from yt.extern.six import string_types
from yt.funcs import \
    validate_width_tuple, \
    fix_unitary, \
    iterable
from yt.units.yt_array import \
    YTArray, YTQuantity
from yt.utilities.exceptions import \
    YTCoordinateNotImplemented, \
    YTInvalidWidthError

def _unknown_coord(field, data):
    raise YTCoordinateNotImplemented

def _get_coord_fields(axi, units = "code_length"):
    def _dds(field, data):
        rv = data.ds.arr(data.fwidth[...,axi].copy(), units)
        return data._reshape_vals(rv)
    def _coords(field, data):
        rv = data.ds.arr(data.fcoords[...,axi].copy(), units)
        return data._reshape_vals(rv)
    return _dds, _coords

def _get_vert_fields(axi, units = "code_length"):
    def _vert(field, data):
        rv = data.ds.arr(data.fcoords_vertex[...,axi].copy(), units)
        return rv
    return _vert

def validate_iterable_width(width, ds, unit=None):
    if isinstance(width[0], tuple) and isinstance(width[1], tuple):
        validate_width_tuple(width[0])
        validate_width_tuple(width[1])
        return (ds.quan(width[0][0], fix_unitary(width[0][1])),
                ds.quan(width[1][0], fix_unitary(width[1][1])))
    elif isinstance(width[0], Number) and isinstance(width[1], Number):
        return (ds.quan(width[0], 'code_length'),
                ds.quan(width[1], 'code_length'))
    elif isinstance(width[0], YTQuantity) and isinstance(width[1], YTQuantity):
        return (ds.quan(width[0]), ds.quan(width[1]))
    else:
        validate_width_tuple(width)
        # If width and unit are both valid width tuples, we
        # assume width controls x and unit controls y
        try:
            validate_width_tuple(unit)
            return (ds.quan(width[0], fix_unitary(width[1])),
                    ds.quan(unit[0], fix_unitary(unit[1])))
        except YTInvalidWidthError:
            return (ds.quan(width[0], fix_unitary(width[1])),
                    ds.quan(width[0], fix_unitary(width[1])))

class CoordinateHandler(object):
    name = None
    
    def __init__(self, ds, ordering):
        self.ds = weakref.proxy(ds)
        self.axis_order = ordering

    def setup_fields(self):
        # This should return field definitions for x, y, z, r, theta, phi
        raise NotImplementedError

    def pixelize(self, dimension, data_source, field, bounds, size, antialias = True):
        # This should *actually* be a pixelize call, not just returning the
        # pixelizer
        raise NotImplementedError

    def distance(self, start, end):
        p1 = self.convert_to_cartesian(start)
        p2 = self.convert_to_cartesian(end)
        return np.sqrt(((p1-p2)**2.0).sum())

    def convert_from_cartesian(self, coord):
        raise NotImplementedError

    def convert_to_cartesian(self, coord):
        raise NotImplementedError

    def convert_to_cylindrical(self, coord):
        raise NotImplementedError

    def convert_from_cylindrical(self, coord):
        raise NotImplementedError

    def convert_to_spherical(self, coord):
        raise NotImplementedError

    def convert_from_spherical(self, coord):
        raise NotImplementedError

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
            rv[i] = (self.axis_name[self.x_axis[i]],
                     self.axis_name[self.y_axis[i]])
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
    def period(self):
        raise NotImplementedError

    def sanitize_depth(self, depth):
        if iterable(depth):
            validate_width_tuple(depth)
            depth = (self.ds.quan(depth[0], fix_unitary(depth[1])), )
        elif isinstance(depth, Number):
            depth = (self.ds.quan(depth, 'code_length',
                                  registry=self.ds.unit_registry), )
        elif isinstance(depth, YTQuantity):
            depth = (depth, )
        else:
            raise YTInvalidWidthError(depth)
        return depth

    def sanitize_width(self, axis, width, depth):
        if width is None:
            # Default to code units
            if not iterable(axis):
                xax = self.x_axis[axis]
                yax = self.y_axis[axis]
                w = self.ds.domain_width[[xax, yax]]
            else:
                # axis is actually the normal vector
                # for an off-axis data object.
                mi = np.argmin(self.ds.domain_width)
                w = self.ds.domain_width[[mi,mi]]
            width = (w[0], w[1])
        elif iterable(width):
            width = validate_iterable_width(width, self.ds)
        elif isinstance(width, YTQuantity):
            width = (width, width)
        elif isinstance(width, Number):
            width = (self.ds.quan(width, 'code_length'),
                     self.ds.quan(width, 'code_length'))
        else:
            raise YTInvalidWidthError(width)
        if depth is not None:
            depth = self.sanitize_depth(depth)
            return width + depth
        return width

    def sanitize_center(self, center, axis):
        if isinstance(center, string_types):
            if center.lower() == "m" or center.lower() == "max":
                v, center = self.ds.find_max(("gas", "density"))
                center = self.ds.arr(center, 'code_length')
            elif center.lower() == "c" or center.lower() == "center":
                center = (self.ds.domain_left_edge + self.ds.domain_right_edge) / 2
            else:
                raise RuntimeError('center keyword \"%s\" not recognized' % center)
        elif isinstance(center, YTArray):
            return self.ds.arr(center), self.convert_to_cartesian(center)
        elif iterable(center):
            if isinstance(center[0], string_types) and isinstance(center[1], string_types):
                if center[0].lower() == "min":
                    v, center = self.ds.find_min(center[1])
                elif center[0].lower() == "max":
                    v, center = self.ds.find_max(center[1])
                else:
                    raise RuntimeError("center keyword \"%s\" not recognized" % center)
                center = self.ds.arr(center, 'code_length')
            elif iterable(center[0]) and isinstance(center[1], string_types):
                center = self.ds.arr(center[0], center[1])
            else:
                center = self.ds.arr(center, 'code_length')
        else:
            raise RuntimeError("center keyword \"%s\" not recognized" % center)
        # This has to return both a center and a display_center
        display_center = self.convert_to_cartesian(center)
        return center, display_center


def cartesian_to_cylindrical(coord, center = (0,0,0)):
    c2 = np.zeros_like(coord)
    if not isinstance(center, YTArray):
        center = center * coord.uq
    c2[...,0] = ((coord[...,0] - center[0])**2.0
              +  (coord[...,1] - center[1])**2.0)**0.5
    c2[...,1] = coord[...,2] # rzt
    c2[...,2] = np.arctan2(coord[...,1] - center[1],
                           coord[...,0] - center[0])
    return c2

def cylindrical_to_cartesian(coord, center = (0,0,0)):
    c2 = np.zeros_like(coord)
    if not isinstance(center, YTArray):
        center = center * coord.uq
    c2[...,0] = np.cos(coord[...,0]) * coord[...,1] + center[0]
    c2[...,1] = np.sin(coord[...,0]) * coord[...,1] + center[1]
    c2[...,2] = coord[...,2]
    return c2

