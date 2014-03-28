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
from yt.utilities.math_utils import get_rotation_matrix
import numpy as np


class Camera(Orientation):

    r"""    """

    _moved = True

    def __init__(self, data_source=None):
        """Initialize a Camera Instance"""
        self.data_source = data_source
        self.position = None
        self.north_vector = None
        self.resolution = (256, 256)
        self.light = None
        self.width = None
        self.focus = np.zeros(3)
        self.position = np.ones(3)
        if data_source is not None:
            self.inherit_default_from_data_source()
        else:
            super(Camera, self).__init__(self.focus - self.position,
                                         self.north_vector, steady_north=False)

    def inherit_default_from_data_source(self):
        data_source = self.data_source
        self.position = data_source.pf.domain_right_edge

        width = 1.5 * data_source.pf.domain_width.max()
        (xmi, xma), (ymi, yma), (zmi, zma) = \
            data_source.quantities['Extrema'](['x', 'y', 'z'])
        width = np.sqrt((xma-xmi)**2 + (yma-ymi)**2 + (zma-zmi)**2) /\
            np.sqrt(3)
        focus = data_source.get_field_parameter('center')

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
                                     self.north_vector, steady_north=False)
        self._moved = True

    def switch_orientation(self, normal_vector=None, north_vector=None):
        r"""Change the view direction based on any of the orientation parameters.

        This will recalculate all the necessary vectors and vector planes related
        to an orientable object.

        Parameters
        ----------
        normal_vector: array_like, optional
            The new looking vector.
        north_vector : array_like, optional
            The 'up' direction for the plane of rays.  If not specific,
            calculated automatically.
        """
        if north_vector is None:
            north_vector = self.north_vector
        if normal_vector is None:
            normal_vector = self.normal_vector
        self._setup_normalized_vectors(normal_vector, north_vector)

    def switch_view(self, normal_vector=None, north_vector=None):
        r"""Change the view based on any of the view parameters.

        This will recalculate the orientation and width based on any of
        normal_vector, width, center, and north_vector.

        Parameters
        ----------
        normal_vector: array_like, optional
            The new looking vector.
        north_vector : array_like, optional
            The 'up' direction for the plane of rays.  If not specific,
            calculated automatically.
        """
        if north_vector is None:
            north_vector = self.north_vector
        if normal_vector is None:
            normal_vector = self.normal_vector
        self.switch_orientation(normal_vector=normal_vector,
                                         north_vector=north_vector)
        self._moved = True

    def pitch(self, theta):
        r"""Rotate by a given angle about the horizontal axis

        Pitch the view.

        Parameters
        ----------
        theta : float, in radians
             Angle (in radians) by which to pitch the view.

        Examples
        --------

        >>> cam.roll(np.pi/4)
        """
        rot_vector = self.unit_vectors[0]
        R = get_rotation_matrix(theta, rot_vector)
        self.switch_view(
            normal_vector=np.dot(R, self.unit_vectors[2]),
            north_vector=np.dot(R, self.unit_vectors[1]))
        if self.steady_north:
            self.north_vector = self.unit_vectors[1]

    def yaw(self, theta):
        r"""Rotate by a given angle about the vertical axis

        Yaw the view.

        Parameters
        ----------
        theta : float, in radians
             Angle (in radians) by which to yaw the view.

        Examples
        --------

        >>> cam.roll(np.pi/4)
        """
        rot_vector = self.unit_vectors[1]
        R = get_rotation_matrix(theta, rot_vector)
        self.switch_view(
            normal_vector=np.dot(R, self.unit_vectors[2]))

    def roll(self, theta):
        r"""Rotate by a given angle about the view normal axis

        Roll the view.

        Parameters
        ----------
        theta : float, in radians
             Angle (in radians) by which to roll the view.

        Examples
        --------

        >>> cam.roll(np.pi/4)
        """
        rot_vector = self.unit_vectors[2]
        R = get_rotation_matrix(theta, rot_vector)
        self.switch_view(
            normal_vector=np.dot(R, self.unit_vectors[2]),
            north_vector=np.dot(R, self.unit_vectors[1]))
        if self.steady_north:
            self.north_vector = np.dot(R, self.north_vector)

