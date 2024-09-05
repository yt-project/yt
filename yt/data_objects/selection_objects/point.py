from yt.data_objects.selection_objects.data_selection_objects import (
    YTSelectionContainer,
    YTSelectionContainer0D,
)
from yt.data_objects.static_output import Dataset
from yt.funcs import validate_3d_array, validate_object
from yt.units import YTArray


class YTPoint(YTSelectionContainer0D):
    """
    A 0-dimensional object defined by a single point

    Parameters
    ----------
    p: array_like
        A points defined within the domain.  If the domain is
        periodic its position will be corrected to lie inside
        the range [DLE,DRE) to ensure one and only one cell may
        match that point
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
    >>> point = ds.point(c)
    """

    _type_name = "point"
    _con_args = ("p",)

    def __init__(self, p, ds=None, field_parameters=None, data_source=None):
        validate_3d_array(p)
        validate_object(ds, Dataset)
        validate_object(field_parameters, dict)
        validate_object(data_source, YTSelectionContainer)
        super().__init__(ds, field_parameters, data_source)
        if isinstance(p, YTArray):
            # we pass p through ds.arr to ensure code units are attached
            self.p = self.ds.arr(p)
        else:
            self.p = self.ds.arr(p, "code_length")
