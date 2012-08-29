"""
Coordinate handler base class.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/
License:
  Copyright (C) 2007-2011 Matthew Turk.  All Rights Reserved.

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
    
    __metaclass__ = abc.ABCMeta

    def __init__(self, pf):
        self.pf = weakref.proxy(pf)

    @abc.abstractmethod
    def coordinate_fields(self):
        # This should return field definitions for x, y, z, r, theta, phi
        pass

    @abc.abstractmethod
    def pixelize(self, dimension, data_source, field, bounds, size, antialias = True):
        # This should *actually* be a pixelize call, not just returning the
        # pixelizer
        pass

    def cartesian_length(self, start, end):
        p1 = self.convert_to_cartesian(start)
        p2 = self.convert_to_cartesian(end)
        return np.sqrt(((p1-p2)**2.0).sum())

    @abc.abstractmethod
    def convert_from_cartesian(self, coord):
        pass

    @abc.abstractmethod
    def convert_to_cartesian(self, coord):
        pass

    @abc.abstractproperty
    def axis_name(self):
        pass

    @abc.abstractproperty
    def axis_id(self):
        pass

    @abc.abstractproperty
    def x_axis(self):
        pass

    @abc.abstractproperty
    def y_axis(self):
        pass

    @abc.abstractproperty
    def period(self):
        pass

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
        buff = _MPL.CPixelize(data_source['x'],   data_source['y'],   data_source['z'],
                              data_source['px'],  data_source['py'],
                              data_source['pdx'], data_source['pdy'], data_source['pdz'],
                              data_source.center, data_source._inv_mat, indices,
                              data_source[item], size[0], size[1], bounds).transpose()
        return buff

    def convert_from_cartesian(self, coord):
        return coord

    def convert_to_cartesian(self, coord):
        return coord

    @property
    def axis_name(self):
        return {
             0  : 'x',  1  : 'y',  2  : 'z',
            'x' : 'x', 'y' : 'y', 'z' : 'z',
            'X' : 'x', 'Y' : 'y', 'Z' : 'z',
        }

    @property
    def axis_id(self):
        return {
            'x' : 0, 'y' : 1, 'z' : 2,
             0  : 0,  1  : 1,  2  : 2,
        }

    @property
    def x_axis(self):
        return {
            'x' : 1, 'y' : 0, 'z' : 0,
             0  : 1,  1  : 0,  2  : 0,
        }

    @property
    def y_axis(self):
        return {
            'x' : 2, 'y' : 2, 'z' : 1,
             0  : 2,  1  : 2,  2  : 1,
        }

    @property
    def period(self):
        return self.pf.domain_width

