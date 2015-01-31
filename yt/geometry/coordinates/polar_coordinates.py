"""
Polar fields




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
from yt.utilities.lib.misc_utilities import \
    pixelize_cylinder

class PolarCoordinateHandler(CoordinateHandler):

    def __init__(self, ds, ordering = 'rtz'):
        if ordering != 'rtz': raise NotImplementedError
        super(PolarCoordinateHandler, self).__init__(ds)

    def setup_fields(self, registry):
        # return the fields for r, z, theta
        registry.add_field("dx", function=_unknown_coord)
        registry.add_field("dy", function=_unknown_coord)
        registry.add_field("x", function=_unknown_coord)
        registry.add_field("y", function=_unknown_coord)

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

        f1, f2 = _get_coord_fields(2) 
        registry.add_field(("index", "dz"), function = f1,
                           display_field = False,
                           units = "code_length")
        registry.add_field(("index", "z"), function = f2,
                           display_field = False,
                           units = "code_length")


        def _CylindricalVolume(field, data):
            return data["dtheta"] * data["r"] * data["dr"] * data["dz"]
        registry.add_field("CellVolume", function=_CylindricalVolume)

    def pixelize(self, dimension, data_source, field, bounds, size, antialias = True):
        ax_name = self.axis_name[dimension]
        if ax_name in ('r', 'theta'):
            return self._ortho_pixelize(data_source, field, bounds, size,
                                        antialias)
        elif ax_name == "z":
            return self._polar_pixelize(data_source, field, bounds, size,
                                        antialias)
        else:
            # Pixelizing along a cylindrical surface is a bit tricky
            raise NotImplementedError


    def _ortho_pixelize(self, data_source, field, bounds, size, antialias):
        buff = _MPL.Pixelize(data_source['px'], data_source['py'],
                             data_source['pdx'], data_source['pdy'],
                             data_source[field], size[0], size[1],
                             bounds, int(antialias),
                             True, self.period).transpose()
        return buff

    def _polar_pixelize(self, data_source, field, bounds, size, antialias):
        # Out bounds here will *always* be what plot window thinks are x0, x1,
        # y0, y1, but which will actually be rmin, rmax, thetamin, thetamax.
        buff = pixelize_cylinder(data_source['r'],
                                 data_source['dr'],
                                 data_source['theta'],
                                 data_source['dtheta'] / 2.0, # half-widths
                                 size, data_source[field], bounds)
        return buff

    axis_name = { 0  : 'r',  1  : 'theta',  2  : 'z',
                 'r' : 'r', 'theta' : 'theta', 'z' : 'z',
                 'R' : 'r', 'Theta' : 'theta', 'Z' : 'z'}

    axis_id = { 'r' : 0, 'theta' : 1, 'z' : 2,
                 0  : 0,  1  : 1,  2  : 2}

    x_axis = { 'r' : 1, 'theta' : 0, 'z' : 0,
                0  : 1,  1  : 0,  2  : 0}

    y_axis = { 'r' : 2, 'theta' : 2, 'z' : 1,
                0  : 2,  1  : 2,  2  : 1}

    _image_axis_name = None
    @property
    def image_axis_name(self):    
        if self._image_axis_name is not None:
            return self._image_axis_name
        # This is the x and y axes labels that get displayed.  For
        # non-Cartesian coordinates, we usually want to override these for
        # Cartesian coordinates, since we transform them.
        rv = {0: ('theta', 'z'),
              1: ('x', 'y'),
              2: ('r', 'z')}
        for i in rv.keys():
            rv[self.axis_name[i]] = rv[i]
            rv[self.axis_name[i].upper()] = rv[i]
        self._image_axis_name = rv
        return rv

    def convert_from_cartesian(self, coord):
        return cartesian_to_cylindrical(coord)

    def convert_to_cartesian(self, coord):
        return cylindrical_to_cartesian(coord)

    def convert_to_cylindrical(self, coord):
        return coord

    def convert_from_cylindrical(self, coord):
        return coord

    def convert_to_spherical(self, coord):
        raise NotImplementedError

    def convert_from_spherical(self, coord):
        raise NotImplementedError

    @property
    def period(self):
        return np.array([0.0, 0.0, 2.0*np.pi])

