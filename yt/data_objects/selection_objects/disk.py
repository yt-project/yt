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


class YTDisk(YTSelectionContainer3D):
    """
    By providing a *center*, a *normal*, a *radius* and a *height* we
    can define a cylinder of any proportion.  Only cells whose centers are
    within the cylinder will be selected.

    Parameters
    ----------
    center : array_like
        coordinate to which the normal, radius, and height all reference
    normal : array_like
        the normal vector defining the direction of lengthwise part of the
        cylinder
    radius : float
        the radius of the cylinder
    height : float
        the distance from the midplane of the cylinder to the top and
        bottom planes
    fields : array of fields, optional
        any fields to be pre-loaded in the cylinder object
    ds: ~yt.data_objects.static_output.Dataset, optional
        An optional dataset to use rather than self.ds
    field_parameters : dictionary
         A dictionary of field parameters than can be accessed by derived
         fields.
    data_source: optional
        Draw the selection from the provided data source rather than
        all data associated with the data_set

    Examples
    --------

    >>> import yt
    >>> ds = yt.load("RedshiftOutput0005")
    >>> c = [0.5, 0.5, 0.5]
    >>> disk = ds.disk(c, [1, 0, 0], (1, "kpc"), (10, "kpc"))
    """

    _type_name = "disk"
    _con_args = ("center", "_norm_vec", "radius", "height")

    def __init__(
        self,
        center,
        normal,
        radius,
        height,
        fields=None,
        ds=None,
        field_parameters=None,
        data_source=None,
    ):
        validate_center(center)
        validate_3d_array(normal)
        validate_float(radius)
        validate_float(height)
        validate_sequence(fields)
        validate_object(ds, Dataset)
        validate_object(field_parameters, dict)
        validate_object(data_source, YTSelectionContainer)
        YTSelectionContainer3D.__init__(self, center, ds, field_parameters, data_source)
        self._norm_vec = np.array(normal) / np.sqrt(np.dot(normal, normal))
        self.set_field_parameter("normal", self._norm_vec)
        self.set_field_parameter("center", self.center)
        self.height = fix_length(height, self.ds)
        self.radius = fix_length(radius, self.ds)
        self._d = -1.0 * np.dot(self._norm_vec, self.center)

    def _get_bbox(self):
        """
        Return the minimum bounding box for the disk.
        """
        # http://www.iquilezles.org/www/articles/diskbbox/diskbbox.htm
        pa = self.center + self._norm_vec * self.height
        pb = self.center - self._norm_vec * self.height
        a = pa - pb
        db = self.radius * np.sqrt(1.0 - a.d * a.d / np.dot(a.d, a.d))
        return np.minimum(pa - db, pb - db), np.maximum(pa + db, pb + db)
