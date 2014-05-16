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

class SpectralCubeCoordinateHandler(CartesianCoordinateHandler):

    def __init__(self, pf):
        super(SpectralCubeCoordinateHandler, self).__init__(pf)

        self.axis_name = {}
        self.axis_id = {}

        for axis, axis_name in zip([pf.lon_axis, pf.lat_axis, pf.spec_axis],
                                   ["Image\ x", "Image\ y", pf.spec_name]):
            lower_ax = "xyz"[axis]
            upper_ax = lower_ax.upper()

            self.axis_name[axis] = axis_name
            self.axis_name[lower_ax] = axis_name
            self.axis_name[upper_ax] = axis_name
            self.axis_name[axis_name] = axis_name

            self.axis_id[lower_ax] = axis
            self.axis_id[axis] = axis
            self.axis_id[axis_name] = axis

        self.default_unit_label = {}
        self.default_unit_label[pf.lon_axis] = "pixel"
        self.default_unit_label[pf.lat_axis] = "pixel"
        self.default_unit_label[pf.spec_axis] = pf.spec_unit

        def _spec_axis(ax, x, y):
            p = (x,y)[ax]
            return [self.pf.pixel2spec(pp).v for pp in p]

        self.axis_field = {}
        self.axis_field[self.pf.spec_axis] = _spec_axis

    def convert_to_cylindrical(self, coord):
        raise NotImplementedError

    def convert_from_cylindrical(self, coord):
        raise NotImplementedError

    x_axis = { 'x' : 1, 'y' : 0, 'z' : 0,
                0  : 1,  1  : 0,  2  : 0}

    y_axis = { 'x' : 2, 'y' : 2, 'z' : 1,
                0  : 2,  1  : 2,  2  : 1}
