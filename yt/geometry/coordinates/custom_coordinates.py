"""
Definitions for nonspatial coordinate systems

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2015, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------


import numpy as np
from .coordinate_handler import \
    CoordinateHandler, \
    _get_coord_fields, \
    _get_vert_fields
import yt.visualization._MPL as _MPL
from collections import OrderedDict


class CustomCoordinateHandler(CoordinateHandler):

    def __init__(self, ds, ordering=(('x', 'code_length'),
                                     ('y', 'code_length'),
                                     ('z', 'code_length'))):
        self.axes_units = OrderedDict(ordering)
        super(CustomCoordinateHandler, self).__init__(ds, tuple(self.axes_units.keys()))

    def setup_fields(self, registry):
        for axi, ax in enumerate(self.axis_order):
            f1, f2 = _get_coord_fields(axi, self.axes_units[ax])
            registry.add_field(("index", "d%s" % ax), function = f1,
                               display_field = False,
                               units = self.axes_units[ax])
            registry.add_field(("index", "path_element_%s" % ax), function = f1,
                               display_field = False,
                               units = self.axes_units[ax])
            registry.add_field(("index", "%s" % ax), function = f2,
                               display_field = False,
                               units = self.axes_units[ax])
            f3 = _get_vert_fields(axi, self.axes_units[ax])
            registry.add_field(("index", "vertex_%s" % ax), function = f3,
                               display_field = False,
                               units = self.axes_units[ax])
        def _cell_volume(field, data):
            rv  = data["index", "d%s" % self.axis_order[0]].copy(order='K')
            rv *= data["index", "d%s" % self.axis_order[1]]
            rv *= data["index", "d%s" % self.axis_order[2]]
            return rv
        registry.add_field(("index", "cell_volume"), function=_cell_volume,
                           display_field=False, units = "*".join(self.axes_units.values()))
        dfl = [("index", "%s" % ax) for ax in self.axis_order]
        dfl += [("index", "d%s" % ax) for ax in self.axis_order]
        dfl += [("index", "cell_volume")]
        registry.check_derived_fields(dfl)

    def pixelize(self, dimension, data_source, field, bounds, size,
                 antialias = True, periodic = True):
        if dimension < 3:
            return self._ortho_pixelize(data_source, field, bounds, size,
                                        antialias, dimension, periodic)
        else:
            return self._oblique_pixelize(data_source, field, bounds, size,
                                          antialias)

    def _ortho_pixelize(self, data_source, field, bounds, size, antialias,
                        dim, periodic):
        # We should be using fcoords
        period = self.period[:2].copy() # dummy here
        period[0] = self.period[self.x_axis[dim]]
        period[1] = self.period[self.y_axis[dim]]
        if hasattr(period, 'in_units'):
            period = period.in_units(list(self.axes_units.values())).d
        buff = _MPL.Pixelize(data_source['px'], data_source['py'],
                             data_source['pdx'], data_source['pdy'],
                             data_source[field], size[0], size[1],
                             bounds, int(antialias),
                             period, int(periodic)).transpose()
        return buff

    def _oblique_pixelize(self, data_source, field, bounds, size, antialias):
        indices = np.argsort(data_source['dx'])[::-1]
        buff = _MPL.CPixelize(data_source['x'], data_source['y'],
                              data_source['z'], data_source['px'],
                              data_source['py'], data_source['pdx'],
                              data_source['pdy'], data_source['pdz'],
                              data_source.center, data_source._inv_mat, indices,
                              data_source[field], size[0], size[1], bounds).transpose()
        return buff

    def convert_from_cartesian(self, coord):
        return coord

    def convert_to_cartesian(self, coord):
        return coord

    def convert_to_cylindrical(self, coord):
        raise NotImplementedError

    def convert_from_cylindrical(self, coord):
        raise NotImplementedError

    def convert_to_spherical(self, coord):
        raise NotImplementedError

    def convert_from_spherical(self, coord):
        raise NotImplementedError

    @property
    def _x_pairs(self):
        return ((self.axis_order[0], self.axis_order[1]),
                (self.axis_order[1], self.axis_order[2]),
                (self.axis_order[2], self.axis_order[0]))

    @property
    def _y_pairs(self):
        return ((self.axis_order[0], self.axis_order[2]),
                (self.axis_order[1], self.axis_order[0]),
                (self.axis_order[2], self.axis_order[1]))

    @property
    def period(self):
        return self.ds.domain_width

