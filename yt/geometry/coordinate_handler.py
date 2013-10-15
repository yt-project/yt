"""
Coordinate handler base class.




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

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
import yt.visualization._MPL as _MPL

def _unknown_coord(field, data):
    raise YTCoordinateNotImplemented

class CoordinateHandler(object):
    
    def __init__(self, pf):
        self.pf = weakref.proxy(pf)

    def setup_fields(self):
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
    c2[...,2] = np.arctan2(coord[...,1] - center[1],
                           coord[...,0] - center[0])
    return c2

def cylindrical_to_cartesian(coord, center = (0,0,0)):
    c2 = np.zeros_like(coord)
    c2[...,0] = np.cos(coord[...,0]) * coord[...,1] + center[0]
    c2[...,1] = np.sin(coord[...,0]) * coord[...,1] + center[1]
    c2[...,2] = coord[...,2]
    return c2

