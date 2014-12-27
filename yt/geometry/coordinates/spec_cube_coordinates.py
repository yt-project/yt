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

    def __init__(self, ds):
        super(SpectralCubeCoordinateHandler, self).__init__(ds)

        self.axis_name = {}
        self.axis_id = {}

        self.default_unit_label = {}
        if ds.lon_name == "X" and ds.lat_name == "Y":
            names = ["x","y"]
        else:
            names = ["Image\ x", "Image\ y"]
            self.default_unit_label[ds.lon_axis] = "pixel"
            self.default_unit_label[ds.lat_axis] = "pixel"
        names.append(ds.spec_name)
        axes = [ds.lon_axis, ds.lat_axis, ds.spec_axis]
        self.default_unit_label[ds.spec_axis] = ds.spec_unit

        for axis, axis_name in zip(axes, names):

            lower_ax = "xyz"[axis]
            upper_ax = lower_ax.upper()

            self.axis_name[axis] = axis_name
            self.axis_name[lower_ax] = axis_name
            self.axis_name[upper_ax] = axis_name
            self.axis_name[axis_name] = axis_name

            self.axis_id[lower_ax] = axis
            self.axis_id[axis] = axis
            self.axis_id[axis_name] = axis

        def _spec_axis(ax, x, y):
            p = (x,y)[ax]
            return [self.ds.pixel2spec(pp).v for pp in p]

        self.axis_field = {}
        self.axis_field[self.ds.spec_axis] = _spec_axis

    def convert_to_cylindrical(self, coord):
        raise NotImplementedError

    def convert_from_cylindrical(self, coord):
        raise NotImplementedError

    x_axis = { 'x' : 1, 'y' : 0, 'z' : 0,
                0  : 1,  1  : 0,  2  : 0}

    y_axis = { 'x' : 2, 'y' : 2, 'z' : 1,
                0  : 2,  1  : 2,  2  : 1}
