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
from yt.data_objects.yt_array import YTArray
from .coordinate_handler import \
    CoordinateHandler, \
    _unknown_coord

#
# Cylindrical fields
#

class CylindricalCoordinateHandler(CoordinateHandler):

    def __init__(self, pf, ordering = 'rzt'):
        if ordering != 'rzt': raise NotImplementedError
        super(CylindricalCoordinateHandler, self).__init__(pf)

    def setup_fields(self, registry):
        # return the fields for r, z, theta
        registry.add_field(("index", "dx"), function=_unknown_coord)
        registry.add_field(("index", "dy"), function=_unknown_coord)
        registry.add_field(("index", "x"), function=_unknown_coord)
        registry.add_field(("index", "y"), function=_unknown_coord)

        def _dr(field, data):
            return np.ones(data.ActiveDimensions, dtype='float64') * data.dds[0]
        registry.add_field(("index", "dr"),
                 function=_dr,
                 display_field=False,
                 validators=[ValidateSpatial(0)])

        def _dz(field, data):
            return np.ones(data.ActiveDimensions, dtype='float64') * data.dds[1]
        registry.add_field(("index", "dz"),
                 function=_dz,
                 display_field=False,
                 validators=[ValidateSpatial(0)])

        def _dtheta(field, data):
            return np.ones(data.ActiveDimensions, dtype='float64') * data.dds[2]
        registry.add_field(("index", "dtheta"),
                 function=_dtheta,
                 display_field=False,
                 validators=[ValidateSpatial(0)])

        def _coordR(field, data):
            dim = data.ActiveDimensions[0]
            return (np.ones(data.ActiveDimensions, dtype='float64')
                           * np.arange(data.ActiveDimensions[0])[:,None,None]
                    +0.5) * data["index", "dr"] + data.LeftEdge[0]
        registry.add_field(("index", "r"),
                 function=_coordR, display_field=False,
                 validators=[ValidateSpatial(0)])

        def _coordZ(field, data):
            dim = data.ActiveDimensions[1]
            return (np.ones(data.ActiveDimensions, dtype='float64')
                           * np.arange(data.ActiveDimensions[1])[None,:,None]
                    +0.5) * data["index", "dz"] + data.LeftEdge[1]
        registry.add_field(("index", "z"),
                 function=_coordZ, display_field=False,
                 validators=[ValidateSpatial(0)])

        def _coordTheta(field, data):
            dim = data.ActiveDimensions[2]
            return (np.ones(data.ActiveDimensions, dtype='float64')
                           * np.arange(data.ActiveDimensions[2])[None,None,:]
                    +0.5) * data["index", "dtheta"] + data.LeftEdge[2]
        registry.add_field(("index", "theta"),
                 function=_coordTheta, display_field=False,
                 validators=[ValidateSpatial(0)])

        def _CylindricalVolume(field, data):
            return data["index", "dtheta"] \
                 * data["index", "r"] \
                 * data["index", "dr"] \
                 * data["index", "dz"]
        registry.add_field(("index", "cell_volume"),
                 function=_CylindricalVolume)


    def pixelize(self, dimension, data_source, field, bounds, size, antialias = True):
        ax_name = self.axis_name[dimension]
        if ax_name in ('r', 'theta'):
            return self._ortho_pixelize(data_source, field, bounds, size,
                                        antialias)
        elif ax_name == "z":
            return self._cyl_pixelize(data_source, field, bounds, size,
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

    def _cyl_pixelize(self, data_source, field, bounds, size, antialias):
        buff = pixelize_cylinder(data_source['r'],
                                 data_source['dr']/2.0,
                                 data_source['theta'],
                                 data_source['dtheta']/2.0,
                                 size[0], data_source[field], bounds[0])
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

