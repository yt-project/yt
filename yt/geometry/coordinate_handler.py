"""
Coordinate handler base class.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt-project.org/
License:
  Copyright (C) 2012 Matthew Turk.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import numpy as np
import abc
import weakref

from yt.funcs import *
from yt.data_objects.field_info_container import \
    NullFunc, FieldInfoContainer
from yt.utilities.io_handler import io_registry
from yt.utilities.logger import ytLogger as mylog
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    ParallelAnalysisInterface, parallel_splitter
from yt.utilities.lib.misc_utilities import \
    pixelize_cylinder
import yt.visualization._MPL

class CoordinatesHandler(object):
    
    def __init__(self, pf):
        self.pf = weakref.proxy(pf)

    def coordinate_fields(self):
        # This should return field definitions for x, y, z, r, theta, phi
        raise NotImplementedError

    def pixelize(self, dimension, data_source, field, bounds, size, antialias = True):
        # This should *actually* be a pixelize call, not just returning the
        # pixelizer
        raise NotImplementedError

    def distance(self, start, end):
        p1 = self.convert_to_cartesian(start)
        p2 = self.convert_to_cartesian(end)
        return np.sqrt(((p1-p2)**2.0).sum())

    def convert_from_cartesian(self, coord):
        raise NotImplementedError

    def convert_to_cartesian(self, coord):
        raise NotImplementedError

    def convert_to_cylindrical(self, coord):
        raise NotImplementedError

    def convert_from_cylindrical(self, coord):
        raise NotImplementedError

    def convert_to_spherical(self, coord):
        raise NotImplementedError

    def convert_from_spherical(self, coord):
        raise NotImplementedError

    @property
    def axis_name(self):
        raise NotImplementedError

    @property
    def axis_id(self):
        raise NotImplementedError

    @property
    def x_axis(self):
        raise NotImplementedError

    @property
    def y_axis(self):
        raise NotImplementedError

    @property
    def period(self):
        raise NotImplementedError

def cartesian_to_cylindrical(coord, center = (0,0,0)):
    c2 = np.zeros_like(coord)
    c2[...,0] = ((coord[...,0] - center[0])**2.0
              +  (coord[...,1] - center[1])**2.0)**0.5
    c2[...,1] = coord[...,2] # rzt
    c2[...,2] = np.arctans(coord[...,1] - center[1],
                           coord[...,0] - center[0])
    return c2

def cylindrical_to_cartesian(coord, center = (0,0,0)):
    c2 = np.zeros_like(coord)
    c2[...,0] = np.cos(coord[...,0]) * coord[...,1] + center[0]
    c2[...,1] = np.sin(coord[...,0]) * coord[...,1] + center[1]
    c2[...,2] = coord[...,2]
    return c2

class CartesianCoordinatesHandler(CoordinatesHandler):

    def __init__(self, pf):
        super(CartesianCoordinatesHandler, self).__init__(pf)


    def coordinate_fields(self):
        pass

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
                              data_source[item], size[0], size[1], bounds).transpose()
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

class PolarCoordinatesHandler(CoordinatesHandler):

    def __init__(self, pf, ordering = 'rzt'):
        if ordering != 'rzt': raise NotImplementedError
        super(PolarCoordinatesHandler, self).__init__(pf)

    def coordinate_fields(self):
        # return the fields for r, z, theta
        pass

    def pixelize(self, dimension):
        raise NotImplementedError
        if dimension in (0, 2):
            return _MPL.Pixelize
        elif dimension == 2:
            return pixelize_cylinder
        elif dimension > 2:
            raise NotImplementedError

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
        return cylindrical_to_cartesian(coord, center)

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
        return na.array([0.0, 0.0, 2.0*np.pi])

