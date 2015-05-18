"""
Volume Rendering Camera Class

"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from yt.funcs import iterable, mylog, ensure_numpy_array
from yt.utilities.orientation import Orientation
from yt.units.yt_array import YTArray
from yt.utilities.math_utils import get_rotation_matrix
from .utils import data_source_or_all
from .lens import lenses
import numpy as np


class Camera(Orientation):

    r"""    """

    _moved = True
    _width = None
    _focus = None
    _position = None
    _resolution = None

    def __init__(self, data_source=None, lens_type='plane-parallel',
                 auto=False):
        """
        Initialize a Camera Instance

        Parameters
        ----------
        data_source: :class:`AMR3DData` or :class:`Dataset`, optional
            This is the source to be rendered, which can be any arbitrary yt
            data object or dataset.
        lens_type: string, optional
            This specifies the type of lens to use for rendering. Current
            options are 'plane-parallel', 'perspective', and 'fisheye'. See
            :class:`yt.visualization.volume_rendering.lens.Lens` for details.
            Default: 'plane-parallel'
        auto: boolean
            If True, build smart defaults using the data source extent. This
            can be time-consuming to iterate over the entire dataset to find
            the positional bounds. Default: False

        Examples
        --------
        >>> cam = Camera(ds)

        """
        self.lens = None
        self.north_vector = None
        self.normal_vector = None
        self.light = None
        self._resolution = (512, 512)
        self._width = 1.0
        self._focus = np.array([0.0]*3)
        self._position = np.array([1.0]*3)
        self.set_lens(lens_type)
        if data_source is not None:
            data_source = data_source_or_all(data_source)
            self._focus = data_source.ds.domain_center
            self._position = data_source.ds.domain_right_edge
            self._width = 1.5*data_source.ds.domain_width
        if auto:
            self.set_defaults_from_data_source(data_source)

        super(Camera, self).__init__(self.focus - self.position,
                                     self.north_vector, steady_north=False)

        # This should be run on-demand if certain attributes are not set.
        self.lens.setup_box_properties(self)

    def position():
        doc = "The position property."

        def fget(self):
            return self._position

        def fset(self, value):
            self._position = value
            self.switch_orientation()

        def fdel(self):
            del self._position
        return locals()
    position = property(**position())

    def width():
        doc = "The width property."

        def fget(self):
            return self._width

        def fset(self, value):
            self._width = value
            self.switch_orientation()

        def fdel(self):
            del self._width
            self._width = None
        return locals()
    width = property(**width())

    def focus():
        doc = "The focus property."

        def fget(self):
            return self._focus

        def fset(self, value):
            self._focus = value
            self.switch_orientation()

        def fdel(self):
            del self._focus
        return locals()
    focus = property(**focus())

    def resolution():
        doc = "The resolution property."

        def fget(self):
            return self._resolution

        def fset(self, value):
            if iterable(value):
                assert (len(value) == 2)
            else:
                value = (value, value)
            self._resolution = value

        def fdel(self):
            del self._resolution
            self._resolution = None
        return locals()
    resolution = property(**resolution())

    def _get_sampler_params(self, render_source):
        lens_params = self.lens._get_sampler_params(self, render_source)
        lens_params.update(width=self.width)
        return lens_params

    def set_lens(self, lens_type):
        if lens_type not in lenses:
            mylog.error("Lens type not available")
            raise RuntimeError()
        self.lens = lenses[lens_type]()
        self.lens.camera = self

    def set_defaults_from_data_source(self, data_source):
        self.position = data_source.pf.domain_right_edge

        width = 1.5 * data_source.pf.domain_width.max()
        (xmi, xma), (ymi, yma), (zmi, zma) = \
            data_source.quantities['Extrema'](['x', 'y', 'z'])
        width = np.sqrt((xma - xmi) ** 2 + (yma - ymi) ** 2 +
                        (zma - zmi) ** 2)
        focus = data_source.get_field_parameter('center')

        if iterable(width) and len(width) > 1 and isinstance(width[1], str):
            width = data_source.pf.quan(width[0], input_units=width[1])
            # Now convert back to code length for subsequent manipulation
            width = width.in_units("code_length")  # .value
        if not iterable(width):
            width = data_source.pf.arr([width, width, width],
                                       input_units='code_length')
            # left/right, top/bottom, front/back
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
        """This must have been created using ds.arr"""
        assert isinstance(width, YTArray), 'Width must be created with ds.arr'
        if isinstance(width, YTArray):
            width = width.in_units('code_length')

        if not iterable(width):
            width = YTArray([width.d]*3, width.units)  # Can't get code units.
        self.width = width
        self.switch_orientation()

    def set_position(self, position, north_vector=None):
          self.position = position
          self.switch_orientation(normal_vector=self.focus - self.position,
                                  north_vector=north_vector)

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
        self.lens.setup_box_properties(self)

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

    def rotate(self, theta, rot_vector=None):
        r"""Rotate by a given angle

        Rotate the view.  If `rot_vector` is None, rotation will occur
        around the `north_vector`.

        Parameters
        ----------
        theta : float, in radians
             Angle (in radians) by which to rotate the view.
        rot_vector  : array_like, optional
            Specify the rotation vector around which rotation will
            occur.  Defaults to None, which sets rotation around
            `north_vector`

        Examples
        --------

        >>> cam.rotate(np.pi/4)
        """
        rotate_all = rot_vector is not None
        if rot_vector is None:
            rot_vector = self.unit_vectors[0]
        else:
            rot_vector = ensure_numpy_array(rot_vector)
            rot_vector = rot_vector/np.linalg.norm(rot_vector)

        R = get_rotation_matrix(theta, rot_vector)

        normal_vector = self.unit_vectors[2]
        normal_vector = normal_vector/np.sqrt((normal_vector**2).sum())

        if rotate_all:
            self.switch_view(
                normal_vector=np.dot(R, normal_vector),
                north_vector=np.dot(R, self.unit_vectors[1]))
        else:
            self.switch_view(normal_vector=np.dot(R, normal_vector))

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
        >>> cam.pitch(np.pi/4)
        """
        self.rotate(theta, rot_vector=self.unit_vectors[0])

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
        >>> cam.yaw(np.pi/4)
        """
        self.rotate(theta, rot_vector=self.unit_vectors[1])

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
        self.rotate(theta, rot_vector=self.unit_vectors[2])

    def rotation(self, theta, n_steps, rot_vector=None):
        r"""Loop over rotate, creating a rotation

        This will rotate `n_steps` until the current view has been
        rotated by an angle `theta`.

        Parameters
        ----------
        theta : float, in radians
            Angle (in radians) by which to rotate the view.
        n_steps : int
            The number of look_at snapshots to make.
        rot_vector  : array_like, optional
            Specify the rotation vector around which rotation will
            occur.  Defaults to None, which sets rotation around the
            original `north_vector`

        Examples
        --------

        >>> for i in cam.rotation(np.pi, 10):
        ...     im = sc.render("rotation_%04i.png" % i)
        """

        dtheta = (1.0*theta)/n_steps
        for i in xrange(n_steps):
            self.rotate(dtheta, rot_vector=rot_vector)
            yield i

    def move_to(self, final, n_steps, exponential=False):
        r"""Loop over a look_at

        This will yield `n_steps` until the current view has been
        moved to a final center of `final`.

        Parameters
        ----------
        final : YTArray
            The final center to move to after `n_steps`
        n_steps : int
            The number of look_at snapshots to make.
        exponential : boolean
            Specifies whether the move/zoom transition follows an
            exponential path toward the destination or linear

        Examples
        --------

        >>> for i in cam.move_to([0.2,0.3,0.6], 10):
        ...     sc.render("move_%04i.png" % i)
        """
        assert isinstance(final, YTArray)
        if exponential:
            position_diff = (np.array(final)/self.position)*1.0
            dx = position_diff**(1.0/n_steps)
        else:
            dx = (final - self.position)*1.0/n_steps
        for i in xrange(n_steps):
            if exponential:
                self.set_position(self.position * dx)
            else:
                self.set_position(self.position + dx)
            yield i

    def zoom(self, factor):
        r"""Change the distance to the focal point.

        This will zoom the camera in by some `factor` toward the focal point,
        along the current view direction, modifying the left/right and up/down
        extents as well.

        Parameters
        ----------
        factor : float
            The factor by which to reduce the distance to the focal point.


        Notes
        -----

        You will need to call snapshot() again to get a new image.

        """
        self.set_width(self.width / factor)

    def zoomin(self, final, n_steps):
        r"""Loop over a zoomin and return snapshots along the way.

        This will yield `n_steps` snapshots until the current view has been
        zooming in to a final factor of `final`.

        Parameters
        ----------
        final : float
            The zoom factor, with respect to current, desired at the end of the
            sequence.
        n_steps : int
            The number of zoom snapshots to make.

        Examples
        --------

        >>> for i in cam.zoomin(100.0, 10):
        ...     sc.render("zoom_%04i.png" % i)
        """
        f = final**(1.0/n_steps)
        for i in xrange(n_steps):
            self.zoom(f)
            yield i

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

    def __repr__(self):
        disp = ("<Camera Object>:\n\tposition:%s\n\tfocus:%s\n\t" +
                "north_vector:%s\n\twidth:%s\n\tlight:%s\n\tresolution:%s\n") \
            % (self.position, self.focus, self.north_vector, self.width,
               self.light, self.resolution)
        disp += "Lens: %s" % self.lens
        return disp
