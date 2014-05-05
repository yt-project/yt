"""
Cartesian fields




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
from .cartesian_coordinates import \
    CartesianCoordinateHandler

class PPVCoordinateHandler(CartesianCoordinateHandler):

    def __init__(self, pf):
        super(PPVCoordinateHandler, self).__init__(pf)

        self.axis_name = {}
        self.axis_id = {}
        self.x_axis = {}
        self.y_axis = {}

        for axis, axis_name in zip([pf.lon_axis, pf.lat_axis, pf.vel_axis],
                                   ["Image\ x", "Image\ y", pf.vel_name]):
            lower_ax = "xyz"[axis]
            upper_ax = lower_ax.upper()

            self.axis_name[axis] = axis_name
            self.axis_name[lower_ax] = axis_name
            self.axis_name[upper_ax] = axis_name
            self.axis_name[axis_name] = axis_name

            self.axis_id[lower_ax] = axis
            self.axis_id[axis] = axis
            self.axis_id[axis_name] = axis

            if axis == 0:
                self.x_axis[axis] = 1
                self.x_axis[lower_ax] = 1
                self.x_axis[axis_name] = 1
            else:
                self.x_axis[axis] = 0
                self.x_axis[lower_ax] = 0
                self.x_axis[axis_name] = 0

            if axis == 2:
                self.y_axis[axis] = 1
                self.y_axis[lower_ax] = 1
                self.y_axis[axis_name] = 1
            else:
                self.y_axis[axis] = 2
                self.y_axis[lower_ax] = 2
                self.y_axis[axis_name] = 2

        self.default_unit_label = {}
        self.default_unit_label[pf.lon_axis] = "pixel"
        self.default_unit_label[pf.lat_axis] = "pixel"
        self.default_unit_label[pf.vel_axis] = pf.vel_unit

    def convert_to_cylindrical(self, coord):
        raise NotImplementedError

    def convert_from_cylindrical(self, coord):
        raise NotImplementedError
