"""
Import the components of the volume rendering extension



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.funcs import iterable
from yt.utilities.orientation import Orientation
from yt.units.yt_array import YTArray
import numpy as np


class Camera(Orientation):

    r"""    """

    _moved = True

    def __init__(self, data_source=None):
        """Initialize a Camera Instance"""
        self.data_source = data_source
        self.position = None
        self.north_vector = np.array([0.0, 0.0, 1.0])
        self.resolution = (256, 256)
        self.light = None
        self.width = None
        self.focus = np.zeros(3)
        self.position = np.ones(3)
        if data_source is not None:
            self.inherit_default_from_data_source()
        else:
            super(Camera, self).__init__(self.focus - self.position,
                                         self.north_vector, steady_north=True)

    def inherit_default_from_data_source(self):
        data_source = self.data_source
        self.position = data_source.pf.domain_right_edge

        width = data_source.pf.domain_width.max()
        focus = data_source.pf.domain_center

        if iterable(width) and len(width) > 1 and isinstance(width[1], str):
            width = self.pf.quan(width[0], input_units=width[1])
            # Now convert back to code length for subsequent manipulation
            width = width.in_units("code_length").value
        if not iterable(width):
            width = (width, width, width)  # left/right, top/bottom, front/back
        if not isinstance(width, YTArray):
            width = self.data_source.pf.arr(width, input_units="code_length")
        if not isinstance(focus, YTArray):
            focus = self.pf.arr(focus, input_units="code_length")

        self.width = width
        self.focus = focus

        super(Camera, self).__init__(self.focus - self.position,
                                     self.north_vector, steady_north=True)
        self._moved = True
