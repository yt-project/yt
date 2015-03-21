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
from .coordinate_handler import \
    _get_coord_fields

class SpectralCubeCoordinateHandler(CartesianCoordinateHandler):

    def __init__(self, ds, ordering = ('x', 'y', 'z')):
        super(SpectralCubeCoordinateHandler, self).__init__(ds, None)

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

        self._image_axis_name = {}

        for axis in axes:

            lower_ax = "xyz"[axis]
            upper_ax = lower_ax.upper()

            ax_names = (names[self.x_axis[axis]], names[self.y_axis[axis]])

            self._image_axis_name[axis] = ax_names
            self._image_axis_name[lower_ax] = ax_names
            self._image_axis_name[upper_ax] = ax_names

        def _spec_axis(ax, x, y):
            p = (x,y)[ax]
            return [self.ds.pixel2spec(pp).v for pp in p]

        self.axis_field = {}
        self.axis_field[self.ds.spec_axis] = _spec_axis

    def setup_fields(self, registry):
        if self.ds.no_cgs_equiv_length == False:
            return super(SpectralCubeCoordinateHandler, self
                    ).setup_fields(registry)
        for axi, ax in enumerate("xyz"):
            f1, f2 = _get_coord_fields(axi)
            def _get_length_func():
                def _length_func(field, data):
                    # Just use axis 0
                    rv = data.ds.arr(data.fcoords[...,0].copy(), field.units)
                    rv[:] = 1.0
                    return rv
                return _length_func
            registry.add_field(("index", "d%s" % ax), function = f1,
                               display_field = False,
                               units = "code_length")
            registry.add_field(("index", "path_element_%s" % ax),
                               function = _get_length_func(),
                               display_field = False,
                               units = "")
            registry.add_field(("index", "%s" % ax), function = f2,
                               display_field = False,
                               units = "code_length")
        def _cell_volume(field, data):
            rv  = data["index", "dx"].copy(order='K')
            rv *= data["index", "dy"]
            rv *= data["index", "dz"]
            return rv
        registry.add_field(("index", "cell_volume"), function=_cell_volume,
                           display_field=False, units = "code_length**3")
        registry.check_derived_fields(
            [("index", "dx"), ("index", "dy"), ("index", "dz"),
             ("index", "x"), ("index", "y"), ("index", "z"),
             ("index", "cell_volume")])

    def convert_to_cylindrical(self, coord):
        raise NotImplementedError

    def convert_from_cylindrical(self, coord):
        raise NotImplementedError

    @property
    def image_axis_name(self):
        return self._image_axis_name

