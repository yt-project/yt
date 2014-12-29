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
from yt.units.yt_array import YTArray
from .coordinate_handler import \
    CoordinateHandler, \
    _unknown_coord, \
    _get_coord_fields
from .cylindrical_coordinates import CylindricalCoordinateHandler
import yt.visualization._MPL as _MPL
from yt.utilities.lib.misc_utilities import \
    pixelize_cylinder

class PolarCoordinateHandler(CylindricalCoordinateHandler):

    def __init__(self, ds, ordering = 'rtz'):
        if ordering != 'rtz': raise NotImplementedError
        super(PolarCoordinateHandler, self).__init__(ds)

    axis_name = { 0  : 'r',  1  : 'theta',  2  : 'z',
                 'r' : 'r', 'theta' : 'theta', 'z' : 'z',
                 'R' : 'r', 'Theta' : 'theta', 'Z' : 'z'}

    axis_id = { 'r' : 0, 'theta' : 1, 'z' : 2,
                 0  : 0,  1  : 1,  2  : 2}

    x_axis = { 'r' : 1, 'theta' : 0, 'z' : 0,
                0  : 1,  1  : 0,  2  : 0}

    y_axis = { 'r' : 2, 'theta' : 2, 'z' : 1,
                0  : 2,  1  : 2,  2  : 1}

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

