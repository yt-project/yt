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
from yt.data_objects.yt_array import YTArray
from .coordinate_handler import \
    CoordinateHandler, \
    _unknown_coord

class CartesianCoordinateHandler(CoordinateHandler):

    def __init__(self, pf):
        super(CartesianCoordinateHandler, self).__init__(pf)

    def setup_fields(self, registry):
        def _get_coord_fields(axi, ax):
            def _dds(field, data):
                return YTArray(data.fwidth[...,axi], 'code_length')
            def _coords(field, data):
                return YTArray(data.fcoords[...,axi], 'code_length')
            return _dds, _coords
        for axi, ax in enumerate('xyz'):
            f1, f2 = _get_coord_fields(axi, ax)
            registry.add_field(("index", "d%s" % ax), function = f1,
                               display_field = False,
                               units = "code_length")
            registry.add_field(("index", "%s" % ax), function = f2,
                               display_field = False,
                               units = "code_length")
        def _cell_volume(field, data):
            return data.fwidth.prod(axis=1)
        registry.add_field(("index", "cell_volume"), function=_cell_volume,
                           display_field=False, units = "code_length**3")
        registry.check_derived_fields(
            [("index", "dx"), ("index", "dy"), ("index", "dz"),
             ("index", "x"), ("index", "y"), ("index", "z"),
             ("index", "cell_volume")])

    def pixelize(self, dimension, data_source, field, bounds, size, antialias = True):
        if dimension < 3:
            return self._ortho_pixelize(data_source, field, bounds, size, antialias)
        else:
            return self._oblique_pixelize(data_source, field, bounds, size, antialias)

    def _ortho_pixelize(self, data_source, field, bounds, size, antialias):
        # We should be using fcoords
        buff = _MPL.Pixelize(data_source['px'], data_source['py'],
                             data_source['pdx'], data_source['pdy'],
                             data_source[field], size[0], size[1],
                             bounds, int(antialias),
                             True, self.period).transpose()
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
        center = self.pf.domain_center
        return cartesian_to_cylindrical(coord, center)

    def convert_from_cylindrical(self, coord):
        center = self.pf.domain_center
        return cylindrical_to_cartesian(coord, center)

    def convert_to_spherical(self, coord):
        raise NotImplementedError

    def convert_from_spherical(self, coord):
        raise NotImplementedError

    # Despite being mutables, we uses these here to be clear about how these
    # are generated and to ensure that they are not re-generated unnecessarily
    axis_name = { 0  : 'x',  1  : 'y',  2  : 'z',
                 'x' : 'x', 'y' : 'y', 'z' : 'z',
                 'X' : 'x', 'Y' : 'y', 'Z' : 'z'}

    axis_id = { 'x' : 0, 'y' : 1, 'z' : 2,
                 0  : 0,  1  : 1,  2  : 2}

    x_axis = { 'x' : 1, 'y' : 0, 'z' : 0,
                0  : 1,  1  : 0,  2  : 0}

    y_axis = { 'x' : 2, 'y' : 2, 'z' : 1,
                0  : 2,  1  : 2,  2  : 1}

    @property
    def period(self):
        return self.pf.domain_width

