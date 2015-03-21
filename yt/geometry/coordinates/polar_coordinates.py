"""
Polar fields




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
from yt.units.yt_array import YTArray
from .coordinate_handler import \
    CoordinateHandler, \
    _unknown_coord, \
    cylindrical_to_cartesian, \
    _get_coord_fields
from .cylindrical_coordinates import CylindricalCoordinateHandler
import yt.visualization._MPL as _MPL
from yt.utilities.lib.misc_utilities import \
    pixelize_cylinder

class PolarCoordinateHandler(CylindricalCoordinateHandler):

  def __init__(self, ds, ordering = ('r', 'theta', 'z')):
        super(PolarCoordinateHandler, self).__init__(ds, ordering)
