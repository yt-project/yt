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
from yt.units.unit_registry import UnitParseError
from yt.utilities.math_utils import get_rotation_matrix
from yt.extern.six import string_types
from .utils import data_source_or_all
from .lens import lenses
import numpy as np


class Camera(Orientation):

    r"""A representation of a point of view into a Scene.

    It is defined by a position (the location of the camera
    in the simulation domain,), a focus (the point at which the
    camera is pointed), a width (the width of the snapshot that will
    be taken, a resolution (the number of pixels in the image), and
    a north_vector (the "up" direction in the resulting image). A
    camera can use a variety of different Lens objects.

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
    
    In this example, the camera is set using defaults that are chosen
    to be reasonable for the argument Dataset.

    >>> import yt
    >>> from yt.visualization.volume_rendering.api import Camera
    >>> ds = yt.load('IsolatedGalaxy/galaxy0030/galaxy0030')
    >>> cam = Camera(ds)

    Here, we set the camera properties manually:

    >>> import yt
    >>> from yt.visualization.volume_rendering.api import Camera
    >>> cam = Camera()
    >>> cam.position = np.array([0.5, 0.5, -1.0])
    >>> cam.focus = np.array([0.5, 0.5, 0.0])
    >>> cam.north_vector = np.array([1.0, 0.0, 0.0])

    Finally, we create a camera with a non-default lens:

    >>> import yt
    >>> from yt.visualization.volume_rendering.api import Camera
    >>> ds = yt.load('IsolatedGalaxy/galaxy0030/galaxy0030')
    >>> cam = Camera(ds, lens_type='perspective')

    """

    _moved = True
    _width = None
    _focus = None
    _position = None
    _resolution = None

    def __init__(self, data_source=None, lens_type='plane-parallel',
                 auto=False):
        """Initialize a Camera Instance"""
        self.lens = None
        self.north_vector = None
        self.normal_vector = None
        self.light = None
        self._resolution = (512, 512)
        self._width = np.array([1.0, 1.0, 1.0])
        self._focus = np.array([0.0]*3)
        self._position = np.array([1.0]*3)
        if data_source is not None:
            data_source = data_source_or_all(data_source)
            self._focus = data_source.ds.domain_center
            self._position = data_source.ds.domain_right_edge
            self._width = data_source.ds.arr(
                [1.5*data_source.ds.domain_width.max()]*3)
            self._domain_center = data_source.ds.domain_center
            self._domain_width = data_source.ds.domain_width
        if auto:
            self.set_defaults_from_data_source(data_source)

        super(Camera, self).__init__(self.focus - self.position,
                                     self.north_vector, steady_north=False)

        self.set_lens(lens_type)

    def position():
        doc = '''The position is the location of the camera in
               the coordinate system of the simulation. This needs
               to be either a YTArray or a numpy array. If it is a 
               numpy array, it is assumed to be in code units. If it
               is a YTArray, it will be converted to code units 
               automatically. '''

        def fget(self):
            return self._position

        def fset(self, value):
            if isinstance(value, YTArray):
                value = value.in_units("code_length")
            self._position = value
            self.normal_vector = self.focus - self._position
            self.switch_orientation()

        def fdel(self):
            del self._position
        return locals()
    position = property(**position())

    def width():
        doc = '''The width of the region that will be seen in the image. 
               This needs to be either a YTArray or a numpy array. If it 
               is a numpy array, it is assumed to be in code units. If it
               is a YTArray, it will be converted to code units automatically. '''

        def fget(self):
            return self._width

        def fset(self, value):
            if isinstance(value, YTArray):
                value = value.in_units("code_length")
            self._width = value
            self.switch_orientation()

        def fdel(self):
            del self._width
            self._width = None
        return locals()
    width = property(**width())

    def focus():
        doc = '''The focus defines the point the Camera is pointed at. This needs
               to be either a YTArray or a numpy array. If it is a 
               numpy array, it is assumed to be in code units. If it
               is a YTArray, it will be converted to code units 
               automatically. '''

        def fget(self):
            return self._focus

        def fset(self, value):
            if isinstance(value, YTArray):
                value = value.in_units("code_length")
            self._focus = value
            self.switch_orientation()

        def fdel(self):
            del self._focus
        return locals()
    focus = property(**focus())

    def resolution():
        doc = '''The resolution is the number of pixels in the image that
               will be produced. '''

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
        r"""Set the lens to be used with this camera.

        Parameters
        ----------

        lens_type : string
            Must be one of the following:
            'plane-parallel'
            'perspective'
            'stereo-perspective'
            'fisheye'
            'spherical'
            'stereo-spherical'

        """
        if lens_type not in lenses:
            mylog.error("Lens type not available")
            raise RuntimeError()
        self.lens = lenses[lens_type]()
        self.lens.set_camera(self)

    def set_defaults_from_data_source(self, data_source):
        """Resets the camera attributes to their default values"""

        position = data_source.ds.domain_right_edge

        width = 1.5 * data_source.ds.domain_width.max()
        (xmi, xma), (ymi, yma), (zmi, zma) = \
            data_source.quantities['Extrema'](['x', 'y', 'z'])
        width = np.sqrt((xma - xmi) ** 2 + (yma - ymi) ** 2 +
                        (zma - zmi) ** 2)
        focus = data_source.get_field_parameter('center')

        if iterable(width) and len(width) > 1 and isinstance(width[1], string_types):
            width = data_source.ds.quan(width[0], input_units=width[1])
            # Now convert back to code length for subsequent manipulation
            width = width.in_units("code_length")  # .value
        if not iterable(width):
            width = data_source.ds.arr([width, width, width],
                                       input_units='code_length')
            # left/right, top/bottom, front/back
        if not isinstance(width, YTArray):
            width = data_source.ds.arr(width, input_units="code_length")
        if not isinstance(focus, YTArray):
            focus = self.ds.arr(focus, input_units="code_length")

        # We can't use the property setters yet, since they rely on attributes
        # that will not be set up until the base class initializer is called.
        # See Issue #1131.
        self._width = width
        self._focus = focus
        self._position = position
        self._domain_center = data_source.ds.domain_center
        self._domain_width = data_source.ds.domain_width

        super(Camera, self).__init__(self.focus - self.position,
                                     self.north_vector, steady_north=False)
        self._moved = True

    def set_width(self, width):
        r"""Set the width of the image that will be produced by this camera.

        Parameters
        ----------

        width : YTQuantity or 3 element YTArray
            The width of the volume rendering in the horizontal, vertical, and
            depth directions. If a scalar, assumes that the width is the same in
            all three directions.
        """
        try:
            width = width.in_units('code_length')
        except (AttributeError, UnitParseError):
            raise ValueError(
                'Volume rendering width must be a YTArray that can be '
                'converted to code units')

        if not iterable(width):
            width = YTArray([width.d]*3, width.units)  # Can't get code units.
        self.width = width
        self.switch_orientation()

    def set_position(self, position, north_vector=None):
        r"""Set the position of the camera.

        Parameters
        ----------

        position : array_like
            The new position
        north_vector : array_like, optional
            The 'up' direction for the plane of rays.  If not specific,
            calculated automatically.

        """

        self.position = position
        self.switch_orientation(normal_vector=self.focus - self.position,
                                north_vector=north_vector)

    def switch_orientation(self, normal_vector=None, north_vector=None):
        r"""Change the view direction based on any of the orientation parameters.

        This will recalculate all the necessary vectors and vector planes
        related to an orientable object.

        Parameters
        ----------
        normal_vector: array_like, optional
            The new looking vector from the camera to the focus.
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
            The new looking vector from the camera to the focus.
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

    def rotate(self, theta, rot_vector=None, rot_center=None):
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
        rot_center  : array_like, optional
            Specifiy the center around which rotation will occur. Defaults
            to None, which sets rotation around the original camera position
            (i.e. the camera position does not change)

        Examples
        --------

        >>> import yt
        >>> import numpy as np
        >>> from yt.visualization.volume_rendering.api import Camera
        >>> cam = Camera()
        >>> # rotate the camera by pi / 4 radians:
        >>> cam.rotate(np.pi/4.0)  
        >>> # rotate the camera about the y-axis instead of cam.north_vector:
        >>> cam.rotate(np.pi/4.0, np.array([0.0, 1.0, 0.0]))  
        >>> # rotate the camera about the origin instead of its own position:
        >>> cam.rotate(np.pi/4.0, rot_center=np.array([0.0, 0.0, 0.0]))  

        """
        rotate_all = rot_vector is not None
        if rot_vector is None:
            rot_vector = self.north_vector
        if rot_center is None:
            rot_center = self._position
        rot_vector = ensure_numpy_array(rot_vector)
        rot_vector = rot_vector/np.linalg.norm(rot_vector)

        new_position = self._position - rot_center
        R = get_rotation_matrix(theta, rot_vector)
        new_position = np.dot(R, new_position) + rot_center

        if (new_position == self._position).all():
            normal_vector = self.unit_vectors[2]
        else:
            normal_vector = rot_center - new_position
        normal_vector = normal_vector/np.sqrt((normal_vector**2).sum())

        if rotate_all:
            self.switch_view(
                normal_vector=np.dot(R, normal_vector),
                north_vector=np.dot(R, self.unit_vectors[1]))
        else:
            self.switch_view(normal_vector=np.dot(R, normal_vector))
        if (new_position != self._position).any(): self.set_position(new_position)

    def pitch(self, theta, rot_center=None):
        r"""Rotate by a given angle about the horizontal axis

        Pitch the view.

        Parameters
        ----------
        theta : float, in radians
             Angle (in radians) by which to pitch the view.
        rot_center  : array_like, optional
            Specifiy the center around which rotation will occur.

        Examples
        --------

        >>> import yt
        >>> import numpy as np
        >>> from yt.visualization.volume_rendering.api import Camera
        >>> cam = Camera()
        >>> # pitch the camera by pi / 4 radians:
        >>> cam.pitch(np.pi/4.0)  
        >>> # pitch the camera about the origin instead of its own position:
        >>> cam.pitch(np.pi/4.0, rot_center=np.array([0.0, 0.0, 0.0]))

        """
        self.rotate(theta, rot_vector=self.unit_vectors[0], rot_center=rot_center)

    def yaw(self, theta, rot_center=None):
        r"""Rotate by a given angle about the vertical axis

        Yaw the view.

        Parameters
        ----------
        theta : float, in radians
             Angle (in radians) by which to yaw the view.
        rot_center  : array_like, optional
            Specifiy the center around which rotation will occur.

        Examples
        --------

        >>> import yt
        >>> import numpy as np
        >>> from yt.visualization.volume_rendering.api import Camera
        >>> cam = Camera()
        >>> # yaw the camera by pi / 4 radians:
        >>> cam.yaw(np.pi/4.0)  
        >>> # yaw the camera about the origin instead of its own position:
        >>> cam.yaw(np.pi/4.0, rot_center=np.array([0.0, 0.0, 0.0]))

        """
        self.rotate(theta, rot_vector=self.unit_vectors[1], rot_center=rot_center)

    def roll(self, theta, rot_center=None):
        r"""Rotate by a given angle about the view normal axis

        Roll the view.

        Parameters
        ----------
        theta : float, in radians
             Angle (in radians) by which to roll the view.
        rot_center  : array_like, optional
            Specifiy the center around which rotation will occur.

        Examples
        --------

        >>> import yt
        >>> import numpy as np
        >>> from yt.visualization.volume_rendering.api import Camera
        >>> cam = Camera()
        >>> # roll the camera by pi / 4 radians:
        >>> cam.roll(np.pi/4.0)  
        >>> # roll the camera about the origin instead of its own position:
        >>> cam.roll(np.pi/4.0, rot_center=np.array([0.0, 0.0, 0.0]))

        """
        self.rotate(theta, rot_vector=self.unit_vectors[2], rot_center=rot_center)

    def iter_rotate(self, theta, n_steps, rot_vector=None, rot_center=None):
        r"""Loop over rotate, creating a rotation

        This will rotate `n_steps` until the current view has been
        rotated by an angle `theta`.

        Parameters
        ----------
        theta : float, in radians
            Angle (in radians) by which to rotate the view.
        n_steps : int
            The number of snapshots to make.
        rot_vector  : array_like, optional
            Specify the rotation vector around which rotation will
            occur.  Defaults to None, which sets rotation around the
            original `north_vector`
        rot_center  : array_like, optional
            Specifiy the center around which rotation will occur. Defaults
            to None, which sets rotation around the original camera position
            (i.e. the camera position does not change)

        Examples
        --------

        >>> import yt
        >>> import numpy as np
        >>> ds = yt.load('IsolatedGalaxy/galaxy0030/galaxy0030')
        >>> 
        >>> im, sc = yt.volume_render(ds)
        >>> cam = sc.camera
        >>> for i in cam.iter_rotate(np.pi, 10):
        ...     im = sc.render()
        ...     sc.save('rotation_%04i.png' % i)

        """

        dtheta = (1.0*theta)/n_steps
        for i in range(n_steps):
            self.rotate(dtheta, rot_vector=rot_vector, rot_center=rot_center)
            yield i

    def iter_move(self, final, n_steps, exponential=False):
        r"""Loop over an iter_move and return snapshots along the way.

        This will yield `n_steps` until the current view has been
        moved to a final center of `final`.

        Parameters
        ----------
        final : YTArray
            The final center to move to after `n_steps`
        n_steps : int
            The number of snapshots to make.
        exponential : boolean
            Specifies whether the move/zoom transition follows an
            exponential path toward the destination or linear.
            Default is False.

        Examples
        --------

        >>> import yt
        >>> import numpy as np
        >>> ds = yt.load('IsolatedGalaxy/galaxy0030/galaxy0030')
        >>> final_position = ds.arr([0.2, 0.3, 0.6], 'unitary')
        >>> im, sc = yt.volume_render(ds)
        >>> cam = sc.camera
        >>> for i in cam.iter_move(final_position, 10):
        ...     sc.render()
        ...     sc.save("move_%04i.png" % i)

        """
        assert isinstance(final, YTArray)
        if exponential:
            position_diff = (final/self.position)*1.0
            dx = position_diff**(1.0/n_steps)
        else:
            dx = (final - self.position)*1.0/n_steps
        for i in range(n_steps):
            if exponential:
                self.set_position(self.position * dx)
            else:
                self.set_position(self.position + dx)
            yield i

    def zoom(self, factor):
        r"""Change the width of the FOV of the camera.

        This will appear to zoom the camera in by some `factor` toward the 
        focal point along the current view direction, but really it's just
        changing the width of the field of view.

        Parameters
        ----------
        factor : float
            The factor by which to divide the width

        Examples
        --------

        >>> import yt
        >>> from yt.visualization.volume_rendering.api import Camera
        >>> ds = yt.load('IsolatedGalaxy/galaxy0030/galaxy0030')
        >>> cam = Camera(ds)
        >>> cam.zoom(1.1)

        """

        self.set_width(self.width / factor)

    def iter_zoom(self, final, n_steps):
        r"""Loop over a iter_zoom and return snapshots along the way.

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

        >>> import yt
        >>> from yt.visualization.volume_rendering.api import Camera
        >>> ds = yt.load('IsolatedGalaxy/galaxy0030/galaxy0030')
        >>> im, sc = yt.volume_render(ds)
        >>> cam = sc.camera
        >>> for i in cam.iter_zoom(100.0, 10):
        ...     sc.render()
        ...     sc.save("zoom_%04i.png" % i)

        """
        f = final**(1.0/n_steps)
        for i in range(n_steps):
            self.zoom(f)
            yield i

    def __repr__(self):
        disp = ("<Camera Object>:\n\tposition:%s\n\tfocus:%s\n\t" +
                "north_vector:%s\n\twidth:%s\n\tlight:%s\n\tresolution:%s\n") \
            % (self.position, self.focus, self.north_vector, self.width,
               self.light, self.resolution)
        disp += "Lens: %s" % self.lens
        return disp
