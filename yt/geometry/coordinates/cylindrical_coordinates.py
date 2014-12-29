"""
Cylindrical fields




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
from yt.units.yt_array import YTArray
from .coordinate_handler import \
    CoordinateHandler, \
    _unknown_coord, \
    _get_coord_fields
import yt.visualization._MPL as _MPL
from yt.utilities.lib.misc_utilities import \
    pixelize_cylinder
#
# Cylindrical fields
#

class CylindricalCoordinateHandler(CoordinateHandler):

    def __init__(self, ds, ordering = 'rzt'):
        if ordering != 'rzt': raise NotImplementedError
        super(CylindricalCoordinateHandler, self).__init__(ds)

    def setup_fields(self, registry):
        # return the fields for r, z, theta
        registry.add_field(("index", "dx"), function=_unknown_coord)
        registry.add_field(("index", "dy"), function=_unknown_coord)
        registry.add_field(("index", "x"), function=_unknown_coord)
        registry.add_field(("index", "y"), function=_unknown_coord)
        f1, f2 = _get_coord_fields(0)
        registry.add_field(("index", "dr"), function = f1,
                           display_field = False,
                           units = "code_length")
        registry.add_field(("index", "r"), function = f2,
                           display_field = False,
                           units = "code_length")

        f1, f2 = _get_coord_fields(self.axis_id['z'])
        registry.add_field(("index", "dz"), function = f1,
                           display_field = False,
                           units = "code_length")
        registry.add_field(("index", "z"), function = f2,
                           display_field = False,
                           units = "code_length")

        f1, f2 = _get_coord_fields(self.axis_id['theta'], "")
        registry.add_field(("index", "dtheta"), function = f1,
                           display_field = False,
                           units = "")
        registry.add_field(("index", "theta"), function = f2,
                           display_field = False,
                           units = "")

        def _CylindricalVolume(field, data):
            return data["index", "dtheta"] \
                 * data["index", "r"] \
                 * data["index", "dr"] \
                 * data["index", "dz"]
        registry.add_field(("index", "cell_volume"),
                 function=_CylindricalVolume,
                 units = "code_length**3")

        def _path_r(field, data):
            return data["index", "dr"]
        registry.add_field(("index", "path_element_r"),
                 function = _path_r,
                 units = "code_length")
        def _path_theta(field, data):
            # Note: this already assumes cell-centered
            return data["index", "r"] * data["index", "dtheta"]
        registry.add_field(("index", "path_element_theta"),
                 function = _path_theta,
                 units = "code_length")
        def _path_z(field, data):
            return data["index", "dz"]
        registry.add_field(("index", "path_element_z"),
                 function = _path_z,
                 units = "code_length")

    def pixelize(self, dimension, data_source, field, bounds, size,
                 antialias = True, periodic = True):
        ax_name = self.axis_name[dimension]
        if ax_name in ('r', 'theta'):
            return self._ortho_pixelize(data_source, field, bounds, size,
                                        antialias, dimension, periodic)
        elif ax_name == "z":
            return self._cyl_pixelize(data_source, field, bounds, size,
                                        antialias)
        else:
            # Pixelizing along a cylindrical surface is a bit tricky
            raise NotImplementedError

    def _ortho_pixelize(self, data_source, field, bounds, size, antialias,
                        dim, periodic):
        period = self.period[:2].copy() # dummy here
        period[0] = self.period[self.x_axis[dim]]
        period[1] = self.period[self.y_axis[dim]]
        if hasattr(period, 'in_units'):
            period = period.in_units("code_length").d
        buff = _MPL.Pixelize(data_source['px'], data_source['py'],
                             data_source['pdx'], data_source['pdy'],
                             data_source[field], size[0], size[1],
                             bounds, int(antialias),
                             period, int(periodic)).transpose()
        return buff

    def _cyl_pixelize(self, data_source, field, bounds, size, antialias):
        buff = pixelize_cylinder(data_source['px'],
                                 data_source['pdx'],
                                 data_source['py'],
                                 data_source['pdy'],
                                 size, data_source[field], bounds)
        return buff

    axis_name = { 0  : 'r',  1  : 'z',  2  : 'theta',
                 'r' : 'r', 'z' : 'z', 'theta' : 'theta',
                 'R' : 'r', 'Z' : 'z', 'Theta' : 'theta'}

    axis_id = { 'r' : 0, 'z' : 1, 'theta' : 2,
                 0  : 0,  1  : 1,  2  : 2}

    x_axis = { 'r' : 1, 'z' : 0, 'theta' : 0,
                0  : 1,  1  : 0,  2  : 0}

    y_axis = { 'r' : 2, 'z' : 2, 'theta' : 1,
                0  : 2,  1  : 2,  2  : 1}

    _image_axis_name = None

    @property
    def image_axis_name(self):    
        if self._image_axis_name is not None:
            return self._image_axis_name
        # This is the x and y axes labels that get displayed.  For
        # non-Cartesian coordinates, we usually want to override these for
        # Cartesian coordinates, since we transform them.
        rv = {0: ('z', 'theta'),
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
