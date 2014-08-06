"""
Spherical fields




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
from .coordinate_handler import \
    CoordinateHandler, \
    _unknown_coord, \
    _get_coord_fields
import yt.visualization._MPL as _MPL
from yt.utilities.lib.misc_utilities import \
    pixelize_cylinder, pixelize_aitoff

class SphericalCoordinateHandler(CoordinateHandler):

    def __init__(self, ds, ordering = 'rtp'):
        if ordering != 'rtp': raise NotImplementedError
        super(SphericalCoordinateHandler, self).__init__(ds)

    def setup_fields(self, registry):
        # return the fields for r, z, theta
        registry.add_field(("index", "dx"), function=_unknown_coord)
        registry.add_field(("index", "dy"), function=_unknown_coord)
        registry.add_field(("index", "dz"), function=_unknown_coord)
        registry.add_field(("index", "x"), function=_unknown_coord)
        registry.add_field(("index", "y"), function=_unknown_coord)
        registry.add_field(("index", "z"), function=_unknown_coord)
        f1, f2 = _get_coord_fields(0)
        registry.add_field(("index", "dr"), function = f1,
                           display_field = False,
                           units = "code_length")
        registry.add_field(("index", "r"), function = f2,
                           display_field = False,
                           units = "code_length")

        f1, f2 = _get_coord_fields(1, "")
        registry.add_field(("index", "dtheta"), function = f1,
                           display_field = False,
                           units = "")
        registry.add_field(("index", "theta"), function = f2,
                           display_field = False,
                           units = "")

        f1, f2 = _get_coord_fields(2, "")
        registry.add_field(("index", "dphi"), function = f1,
                           display_field = False,
                           units = "")
        registry.add_field(("index", "phi"), function = f2,
                           display_field = False,
                           units = "")

        def _SphericalVolume(field, data):
            # r**2 sin theta dr dtheta dphi
            vol = data["index", "r"]**2.0
            vol *= data["index", "dr"]
            vol *= np.sin(data["index", "theta"])
            vol *= data["index", "dtheta"]
            vol *= data["index", "dphi"]
            return vol
        registry.add_field(("index", "cell_volume"),
                 function=_SphericalVolume,
                 units = "code_length**3")

    def pixelize(self, dimension, data_source, field, bounds, size,
                 antialias = True, periodic = True):
        self.period
        if dimension == 0:
            return self._ortho_pixelize(data_source, field, bounds, size,
                                        antialias, dimension, periodic)
        elif dimension in (1, 2):
            return self._cyl_pixelize(data_source, field, bounds, size,
                                          antialias, dimension)
        else:
            raise NotImplementedError

    def _ortho_pixelize(self, data_source, field, bounds, size, antialias,
                        dim, periodic):
        # We should be using fcoords
        period = self.period[:2].copy() # dummy here
        period[0] = self.period[self.x_axis[dim]]
        period[1] = self.period[self.y_axis[dim]]
        period = period.in_units("code_length").d
        buff = _MPL.Pixelize(data_source['px'], data_source['py'],
                             data_source['pdx'], data_source['pdy'],
                             data_source[field], size[0], size[1],
                             bounds, int(antialias),
                             period, int(periodic)).transpose()
        return buff

    def _cyl_pixelize(self, data_source, field, bounds, size, antialias,
                      dimension):
        if dimension == 1:
            buff = pixelize_cylinder(data_source['r'],
                                     data_source['dr'] / 2.0,
                                     data_source['phi'],
                                     data_source['dphi'] / 2.0, # half-widths
                                     size, data_source[field], bounds)
        elif dimension == 2:
            buff = pixelize_cylinder(data_source['r'],
                                     data_source['dr'] / 2.0,
                                     data_source['theta'],
                                     data_source['dtheta'] / 2.0, # half-widths
                                     size, data_source[field], bounds)
            buff = buff.transpose()
        else:
            raise RuntimeError
        return buff


    def convert_from_cartesian(self, coord):
        raise NotImplementedError

    def convert_to_cartesian(self, coord):
        if isinstance(coord, np.ndarray) and len(coord.shape) > 1:
            r = coord[:,0]
            theta = coord[:,1]
            phi = coord[:,2]
            nc = np.zeros_like(coord)
            # r, theta, phi
            nc[:,0] = np.cos(phi) * np.sin(theta)*r
            nc[:,1] = np.sin(phi) * np.sin(theta)*r
            nc[:,2] = np.cos(theta) * r
        else:
            r, theta, phi = coord
            nc = (np.cos(phi) * np.sin(theta)*r,
                  np.sin(phi) * np.sin(theta)*r,
                  np.cos(theta) * r)
        return nc

    def convert_to_cylindrical(self, coord):
        raise NotImplementedError

    def convert_from_cylindrical(self, coord):
        raise NotImplementedError

    def convert_to_spherical(self, coord):
        raise NotImplementedError

    def convert_from_spherical(self, coord):
        raise NotImplementedError

    # Despite being mutables, we uses these here to be clear about how these
    # are generated and to ensure that they are not re-generated unnecessarily
    axis_name = { 0  : 'r',  1  : 'theta',  2  : 'phi',
                 'r' : 'r', 'theta' : 'theta', 'phi' : 'phi',
                 'R' : 'r', 'Theta' : 'theta', 'Phi' : 'phi'}

    _image_axis_name = None
    @property
    def image_axis_name(self):    
        if self._image_axis_name is not None:
            return self._image_axis_name
        # This is the x and y axes labels that get displayed.  For
        # non-Cartesian coordinates, we usually want to override these for
        # Cartesian coordinates, since we transform them.
        rv = {0: ('theta', 'phi'),
              1: ('x / \\sin(\\theta)', 'y / \\sin(\\theta)'),
              2: ('R', 'z')}
        for i in rv.keys():
            rv[self.axis_name[i]] = rv[i]
            rv[self.axis_name[i].upper()] = rv[i]
        self._image_axis_name = rv
        return rv


    axis_id = { 'r' : 0, 'theta' : 1, 'phi' : 2,
                 0  : 0,  1  : 1,  2  : 2}

    x_axis = { 'r' : 1, 'theta' : 0, 'phi' : 0,
                0  : 1,  1  : 0,  2  : 0}

    y_axis = { 'r' : 2, 'theta' : 2, 'phi' : 1,
                0  : 2,  1  : 2,  2  : 1}

    @property
    def period(self):
        return self.ds.domain_width

    def sanitize_center(self, center, axis):
        center, display_center = super(
            SphericalCoordinateHandler, self).sanitize_center(center, axis)
        if axis == 0:
            display_center = center
        elif axis == 1:
            display_center = (0.0 * display_center[0],
                              0.0 * display_center[1],
                              0.0 * display_center[2])
        elif axis ==2:
            display_center = (self.ds.domain_width[0]/2.0,
                              0.0 * display_center[1],
                              0.0 * display_center[2])
        return center, display_center

    def sanitize_width(self, axis, width, depth):
        if width is not None:
            width = super(SphericalCoordinateHandler, self).sanitize_width(
              axis, width, depth)
        elif axis == 0:
            width = [self.ds.domain_width[self.x_axis[0]],
                     self.ds.domain_width[self.y_axis[0]]]
        elif axis == 1:
            # Remember, in spherical coordinates when we cut in theta,
            # we create a conic section
            width = [2.0*self.ds.domain_width[0],
                     2.0*self.ds.domain_width[0]]
        elif axis == 2:
            width = [self.ds.domain_width[0],
                     2.0*self.ds.domain_width[0]]
        return width
