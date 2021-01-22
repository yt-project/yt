from yt.data_objects.selection_objects.data_selection_objects import (
    YTSelectionContainer,
    YTSelectionContainer3D,
)
from yt.data_objects.static_output import Dataset
from yt.funcs import (
    validate_3d_array,
    validate_center,
    validate_object,
    validate_sequence,
)
from yt.units import YTArray


class YTRegion(YTSelectionContainer3D):
    """A 3D region of data with an arbitrary center.

    Takes an array of three *left_edge* coordinates, three
    *right_edge* coordinates, and a *center* that can be anywhere
    in the domain. If the selected region extends past the edges
    of the domain, no data will be found there, though the
    object's `left_edge` or `right_edge` are not modified.

    Parameters
    ----------
    center : array_like
        The center of the region
    left_edge : array_like
        The left edge of the region
    right_edge : array_like
        The right edge of the region
    """

    _type_name = "region"
    _con_args = ("center", "left_edge", "right_edge")

    def __init__(
        self,
        center,
        left_edge,
        right_edge,
        fields=None,
        ds=None,
        field_parameters=None,
        data_source=None,
    ):
        if center is not None:
            validate_center(center)
        validate_3d_array(left_edge)
        validate_3d_array(right_edge)
        validate_sequence(fields)
        validate_object(ds, Dataset)
        validate_object(field_parameters, dict)
        validate_object(data_source, YTSelectionContainer)
        YTSelectionContainer3D.__init__(self, center, ds, field_parameters, data_source)
        if not isinstance(left_edge, YTArray):
            self.left_edge = self.ds.arr(left_edge, "code_length", dtype="float64")
        else:
            # need to assign this dataset's unit registry to the YTArray
            self.left_edge = self.ds.arr(left_edge.copy(), dtype="float64")
        if not isinstance(right_edge, YTArray):
            self.right_edge = self.ds.arr(right_edge, "code_length", dtype="float64")
        else:
            # need to assign this dataset's unit registry to the YTArray
            self.right_edge = self.ds.arr(right_edge.copy(), dtype="float64")

    def _get_bbox(self):
        """
        Return the minimum bounding box for the region.
        """
        return self.left_edge.copy(), self.right_edge.copy()
