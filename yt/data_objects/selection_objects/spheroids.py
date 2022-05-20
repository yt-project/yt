import numpy as np

from yt.data_objects.selection_objects.data_selection_objects import (
    YTSelectionContainer,
    YTSelectionContainer3D,
)
from yt.data_objects.static_output import Dataset
from yt.funcs import (
    fix_length,
    validate_3d_array,
    validate_center,
    validate_float,
    validate_object,
    validate_sequence,
)
from yt.units import YTArray
from yt.utilities.exceptions import YTEllipsoidOrdering, YTException, YTSphereTooSmall
from yt.utilities.logger import ytLogger as mylog
from yt.utilities.math_utils import get_rotation_matrix
from yt.utilities.on_demand_imports import _miniball


class YTSphere(YTSelectionContainer3D):
    """
    A sphere of points defined by a *center* and a *radius*.

    Parameters
    ----------
    center : array_like
        The center of the sphere.
    radius : float, width specifier, or YTQuantity
        The radius of the sphere. If passed a float,
        that will be interpreted in code units. Also
        accepts a (radius, unit) tuple or YTQuantity
        instance with units attached.

    Examples
    --------

    >>> import yt
    >>> ds = yt.load("RedshiftOutput0005")
    >>> c = [0.5, 0.5, 0.5]
    >>> sphere = ds.sphere(c, (1.0, "kpc"))
    """

    _type_name = "sphere"
    _con_args = ("center", "radius")

    def __init__(
        self, center, radius, ds=None, field_parameters=None, data_source=None
    ):
        validate_center(center)
        validate_float(radius)
        validate_object(ds, Dataset)
        validate_object(field_parameters, dict)
        validate_object(data_source, YTSelectionContainer)
        super().__init__(center, ds, field_parameters, data_source)
        # Unpack the radius, if necessary
        radius = fix_length(radius, self.ds)
        if radius < self.index.get_smallest_dx():
            raise YTSphereTooSmall(
                ds,
                radius.in_units("code_length"),
                self.index.get_smallest_dx().in_units("code_length"),
            )
        self.set_field_parameter("radius", radius)
        self.set_field_parameter("center", self.center)
        self.radius = radius

    def _get_bbox(self):
        """
        Return the minimum bounding box for the sphere.
        """
        return -self.radius + self.center, self.radius + self.center


class YTMinimalSphere(YTSelectionContainer3D):
    """
    Build the smallest sphere that encompasses a set of points.

    Parameters
    ----------
    points : YTArray
        The points that the sphere will contain.

    Examples
    --------

    >>> import yt
    >>> ds = yt.load("output_00080/info_00080.txt")
    >>> points = ds.r["particle_position"]
    >>> sphere = ds.minimal_sphere(points)
    """

    _type_name = "sphere"
    _override_selector_name = "minimal_sphere"
    _con_args = ("center", "radius")

    def __init__(self, points, ds=None, field_parameters=None, data_source=None):
        validate_object(ds, Dataset)
        validate_object(field_parameters, dict)
        validate_object(data_source, YTSelectionContainer)
        validate_object(points, YTArray)

        points = fix_length(points, ds)
        if len(points) < 2:
            raise YTException(
                f"Not enough points. Expected at least 2, got {len(points)}"
            )
        mylog.debug("Building minimal sphere around points.")
        mb = _miniball.Miniball(points)
        if not mb.is_valid():
            raise YTException("Could not build valid sphere around points.")

        center = ds.arr(mb.center(), points.units)
        radius = ds.quan(np.sqrt(mb.squared_radius()), points.units)
        super().__init__(center, ds, field_parameters, data_source)
        self.set_field_parameter("radius", radius)
        self.set_field_parameter("center", self.center)
        self.radius = radius


class YTEllipsoid(YTSelectionContainer3D):
    """
    By providing a *center*,*A*,*B*,*C*,*e0*,*tilt* we
    can define a ellipsoid of any proportion.  Only cells whose
    centers are within the ellipsoid will be selected.

    Parameters
    ----------
    center : array_like
        The center of the ellipsoid.
    A : float
        The magnitude of the largest axis (semi-major) of the ellipsoid.
    B : float
        The magnitude of the medium axis (semi-medium) of the ellipsoid.
    C : float
        The magnitude of the smallest axis (semi-minor) of the ellipsoid.
    e0 : array_like (automatically normalized)
        the direction of the largest semi-major axis of the ellipsoid
    tilt : float
        After the rotation about the z-axis to align e0 to x in the x-y
        plane, and then rotating about the y-axis to align e0 completely
        to the x-axis, tilt is the angle in radians remaining to
        rotate about the x-axis to align both e1 to the y-axis and e2 to
        the z-axis.
    Examples
    --------

    >>> import yt
    >>> ds = yt.load("RedshiftOutput0005")
    >>> c = [0.5, 0.5, 0.5]
    >>> ell = ds.ellipsoid(c, 0.1, 0.1, 0.1, np.array([0.1, 0.1, 0.1]), 0.2)
    """

    _type_name = "ellipsoid"
    _con_args = ("center", "_A", "_B", "_C", "_e0", "_tilt")

    def __init__(
        self,
        center,
        A,
        B,
        C,
        e0,
        tilt,
        fields=None,
        ds=None,
        field_parameters=None,
        data_source=None,
    ):
        validate_center(center)
        validate_float(A)
        validate_float(B)
        validate_float(C)
        validate_3d_array(e0)
        validate_float(tilt)
        validate_sequence(fields)
        validate_object(ds, Dataset)
        validate_object(field_parameters, dict)
        validate_object(data_source, YTSelectionContainer)
        YTSelectionContainer3D.__init__(self, center, ds, field_parameters, data_source)
        # make sure the magnitudes of semi-major axes are in order
        if A < B or B < C:
            raise YTEllipsoidOrdering(ds, A, B, C)
        # make sure the smallest side is not smaller than dx
        self._A = self.ds.quan(A, "code_length")
        self._B = self.ds.quan(B, "code_length")
        self._C = self.ds.quan(C, "code_length")
        if self._C < self.index.get_smallest_dx():
            raise YTSphereTooSmall(self.ds, self._C, self.index.get_smallest_dx())
        self._e0 = e0 = e0 / (e0**2.0).sum() ** 0.5
        self._tilt = tilt

        # find the t1 angle needed to rotate about z axis to align e0 to x
        t1 = np.arctan(e0[1] / e0[0])
        # rotate e0 by -t1
        RZ = get_rotation_matrix(t1, (0, 0, 1)).transpose()
        r1 = (e0 * RZ).sum(axis=1)
        # find the t2 angle needed to rotate about y axis to align e0 to x
        t2 = np.arctan(-r1[2] / r1[0])
        """
        calculate the original e1
        given the tilt about the x axis when e0 was aligned
        to x after t1, t2 rotations about z, y
        """
        RX = get_rotation_matrix(-tilt, (1, 0, 0)).transpose()
        RY = get_rotation_matrix(-t2, (0, 1, 0)).transpose()
        RZ = get_rotation_matrix(-t1, (0, 0, 1)).transpose()
        e1 = ((0, 1, 0) * RX).sum(axis=1)
        e1 = (e1 * RY).sum(axis=1)
        e1 = (e1 * RZ).sum(axis=1)
        e2 = np.cross(e0, e1)

        self._e1 = e1
        self._e2 = e2

        self.set_field_parameter("A", A)
        self.set_field_parameter("B", B)
        self.set_field_parameter("C", C)
        self.set_field_parameter("e0", e0)
        self.set_field_parameter("e1", e1)
        self.set_field_parameter("e2", e2)

    def _get_bbox(self):
        """
        Get the bounding box for the ellipsoid. NOTE that in this case
        it is not the *minimum* bounding box.
        """
        radius = self.ds.arr(np.max([self._A, self._B, self._C]), "code_length")
        return -radius + self.center, radius + self.center
