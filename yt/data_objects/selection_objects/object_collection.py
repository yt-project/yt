import numpy as np

from yt.data_objects.selection_objects.data_selection_objects import (
    YTSelectionContainer,
    YTSelectionContainer3D,
)
from yt.data_objects.static_output import Dataset
from yt.funcs import validate_center, validate_object, validate_sequence


class YTDataCollection(YTSelectionContainer3D):
    """
    By selecting an arbitrary *object_list*, we can act on those grids.
    Child cells are not returned.
    """

    _type_name = "data_collection"
    _con_args = ("_obj_list",)

    def __init__(
        self, obj_list, ds=None, field_parameters=None, data_source=None, center=None
    ):
        validate_sequence(obj_list)
        validate_object(ds, Dataset)
        validate_object(field_parameters, dict)
        validate_object(data_source, YTSelectionContainer)
        if center is not None:
            validate_center(center)
        YTSelectionContainer3D.__init__(self, center, ds, field_parameters, data_source)
        self._obj_ids = np.array([o.id - o._id_offset for o in obj_list], dtype="int64")
        self._obj_list = obj_list
