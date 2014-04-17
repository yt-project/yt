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

from yt.funcs import iterable, mylog
from yt.utilities.orientation import Orientation
from yt.units.yt_array import YTArray
from yt.utilities.math_utils import get_rotation_matrix
from utils import data_source_or_all
from lens import lenses
import numpy as np


class Camera(Orientation):

    r"""    """

    _moved = True

    def __init__(self, data_source=None, lens_type='plane-parallel'):
        mylog.debug("Entering %s" % str(self))
        """Initialize a Camera Instance"""
        self.lens = None
        self.position = None
        self.north_vector = None
        self.resolution = (256, 256)
        self.light = None
        self.width = None
        self.focus = np.zeros(3)
        self.position = np.ones(3)
        self.set_lens(lens_type)
        if data_source is not None:
            data_source = data_source_or_all(data_source)
            self.set_defaults_from_data_source(data_source)

        super(Camera, self).__init__(self.focus - self.position,
                                     self.north_vector, steady_north=False)

    def get_sampler_params(self, render_source):
        lens_params = self.lens.get_sampler_params(self, render_source)
        lens_params.update(width=self.width)
        return lens_params

    def set_lens(self, lens_type):
        if lens_type not in lenses:
            mylog.error("Lens type not available")
            raise RuntimeError()
        self.lens = lenses[lens_type]()

    def set_defaults_from_data_source(self, data_source):
        self.position = data_source.pf.domain_right_edge

        width = 1.5 * data_source.pf.domain_width.max()
        (xmi, xma), (ymi, yma), (zmi, zma) = \
            data_source.quantities['Extrema'](['x', 'y', 'z'])
        width = np.sqrt((xma - xmi) ** 2 + (yma - ymi) ** 2 +
                        (zma - zmi) ** 2) / np.sqrt(3)
        focus = data_source.get_field_parameter('center')

        if iterable(width) and len(width) > 1 and isinstance(width[1], str):
            width = data_source.pf.quan(width[0], input_units=width[1])
            # Now convert back to code length for subsequent manipulation
            width = width.in_units("code_length").value
        if not iterable(width):
            width = (width, width, width)  # left/right, top/bottom, front/back
        if not isinstance(width, YTArray):
            width = data_source.pf.arr(width, input_units="code_length")
        if not isinstance(focus, YTArray):
            focus = self.pf.arr(focus, input_units="code_length")

        self.width = width
        self.focus = focus

        super(Camera, self).__init__(self.focus - self.position,
                                     self.north_vector, steady_north=False)
        self._moved = True

    def set_width(self, width):
        if not iterable(width):
            width = [width, width, width] # No way to get code units.
        self.width = width
        self.switch_orientation()

    def set_position(self, position):
        self.position = position
        self.switch_orientation()

    def switch_orientation(self, normal_vector=None, north_vector=None):
        r"""
        Change the view direction based on any of the orientation parameters.

        This will recalculate all the necessary vectors and vector planes
        related to an orientable object.

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

        >>> cam = Camera()
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

        >>> cam = Camera()
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

        >>> cam = Camera()
        >>> cam.roll(np.pi/4)
        """
        rot_vector = self.unit_vectors[2]
        R = get_rotation_matrix(theta, rot_vector)
        self.switch_view(
            normal_vector=np.dot(R, self.unit_vectors[2]),
            north_vector=np.dot(R, self.unit_vectors[1]))
        if self.steady_north:
            self.north_vector = np.dot(R, self.north_vector)

    def project_to_plane(self, pos, res=None):
        if res is None:
            res = self.resolution
        dx = np.dot(pos - self.position.d, self.unit_vectors[1])
        dy = np.dot(pos - self.position.d, self.unit_vectors[0])
        dz = np.dot(pos - self.position.d, self.unit_vectors[2])
        # Transpose into image coords.
        py = (res[0]/2 + res[0]*(dx/self.width[0].d)).astype('int')
        px = (res[1]/2 + res[1]*(dy/self.width[1].d)).astype('int')
        return px, py, dz
