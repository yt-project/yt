"""
Definitions for polar coordinate systems




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from .cylindrical_coordinates import CylindricalCoordinateHandler


class PolarCoordinateHandler(CylindricalCoordinateHandler):
    name = "polar"

    def __init__(self, ds, ordering = ('r', 'theta', 'z')):
        super(PolarCoordinateHandler, self).__init__(ds, ordering)
        # No need to set labels here
