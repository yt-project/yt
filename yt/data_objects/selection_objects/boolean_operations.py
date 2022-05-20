import numpy as np
from more_itertools import always_iterable

import yt.geometry
from yt.data_objects.selection_objects.data_selection_objects import (
    YTSelectionContainer,
    YTSelectionContainer3D,
)
from yt.data_objects.static_output import Dataset
from yt.funcs import validate_object, validate_sequence


class YTBooleanContainer(YTSelectionContainer3D):
    r"""
    This is a boolean operation, accepting AND, OR, XOR, and NOT for combining
    multiple data objects.

    This object is not designed to be created directly; it is designed to be
    created implicitly by using one of the bitwise operations (&, \|, ^, \~) on
    one or two other data objects.  These correspond to the appropriate boolean
    operations, and the resultant object can be nested.

    Parameters
    ----------
    op : string
        Can be AND, OR, XOR, NOT or NEG.
    dobj1 : yt.data_objects.selection_objects.data_selection_objects.
            YTSelectionContainer
        The first selection object
    dobj2 : yt.data_objects.selection_objects.base_objects.YTSelectionContainer
        The second object

    Examples
    --------

    >>> import yt
    >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
    >>> sp = ds.sphere("c", 0.1)
    >>> dd = ds.r[:, :, :]
    >>> new_obj = sp ^ dd
    >>> print(new_obj.sum("cell_volume"), dd.sum("cell_volume") - sp.sum("cell_volume"))
    """

    _type_name = "bool"
    _con_args = ("op", "dobj1", "dobj2")

    def __init__(
        self, op, dobj1, dobj2, ds=None, field_parameters=None, data_source=None
    ):
        YTSelectionContainer3D.__init__(self, None, ds, field_parameters, data_source)
        self.op = op.upper()
        self.dobj1 = dobj1
        self.dobj2 = dobj2
        name = f"Boolean{self.op}Selector"
        sel_cls = getattr(yt.geometry.selection_routines, name)
        self._selector = sel_cls(self)

    def _get_bbox(self):
        le1, re1 = self.dobj1._get_bbox()
        if self.op == "NOT":
            return le1, re1
        else:
            le2, re2 = self.dobj2._get_bbox()
            return np.minimum(le1, le2), np.maximum(re1, re2)


class YTIntersectionContainer3D(YTSelectionContainer3D):
    """
    This is a more efficient method of selecting the intersection of multiple
    data selection objects.

    Creating one of these objects returns the intersection of all of the
    sub-objects; it is designed to be a faster method than chaining & ("and")
    operations to create a single, large intersection.

    Parameters
    ----------
    data_objects : Iterable of YTSelectionContainer
        The data objects to intersect

    Examples
    --------

    >>> import yt
    >>> ds = yt.load("RedshiftOutput0005")
    >>> sp1 = ds.sphere((0.4, 0.5, 0.6), 0.15)
    >>> sp2 = ds.sphere((0.38, 0.51, 0.55), 0.1)
    >>> sp3 = ds.sphere((0.35, 0.5, 0.6), 0.15)
    >>> new_obj = ds.intersection((sp1, sp2, sp3))
    >>> print(new_obj.sum("cell_volume"))
    """

    _type_name = "intersection"
    _con_args = ("data_objects",)

    def __init__(self, data_objects, ds=None, field_parameters=None, data_source=None):
        validate_sequence(data_objects)
        for obj in data_objects:
            validate_object(obj, YTSelectionContainer)
        validate_object(ds, Dataset)
        validate_object(field_parameters, dict)
        validate_object(data_source, YTSelectionContainer)
        YTSelectionContainer3D.__init__(self, None, ds, field_parameters, data_source)
        self.data_objects = list(always_iterable(data_objects))


class YTDataObjectUnion(YTSelectionContainer3D):
    """
    This is a more efficient method of selecting the union of multiple
    data selection objects.

    Creating one of these objects returns the union of all of the sub-objects;
    it is designed to be a faster method than chaining | (or) operations to
    create a single, large union.

    Parameters
    ----------
    data_objects : Iterable of YTSelectionContainer
        The data objects to union

    Examples
    --------

    >>> import yt
    >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
    >>> sp1 = ds.sphere((0.4, 0.5, 0.6), 0.1)
    >>> sp2 = ds.sphere((0.3, 0.5, 0.15), 0.1)
    >>> sp3 = ds.sphere((0.5, 0.5, 0.9), 0.1)
    >>> new_obj = ds.union((sp1, sp2, sp3))
    >>> print(new_obj.sum("cell_volume"))
    """

    _type_name = "union"
    _con_args = ("data_objects",)

    def __init__(self, data_objects, ds=None, field_parameters=None, data_source=None):
        validate_sequence(data_objects)
        for obj in data_objects:
            validate_object(obj, YTSelectionContainer)
        validate_object(ds, Dataset)
        validate_object(field_parameters, dict)
        validate_object(data_source, YTSelectionContainer)
        YTSelectionContainer3D.__init__(self, None, ds, field_parameters, data_source)
        self.data_objects = list(always_iterable(data_objects))
