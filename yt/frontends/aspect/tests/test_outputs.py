from yt.testing import assert_array_equal, assert_equal
from yt.utilities.answer_testing.framework import (
    AllFieldValuesTest,
    can_run_ds,
    data_dir_load,
    requires_ds,
)

_node_names = ("velocity", "p", "T")
_fields = ("T", "p", "velocity_x", "velocity_z", "velocity_y")
box_convect = "aspect/output_convection_box_3d/nproc_1/solution/solution-00002.pvtu"
box_convect_par = "aspect/output_convection_box_3d/nproc_4/solution/solution-00001.pvtu"


@requires_ds(box_convect)
def test_box_convect():
    ds = data_dir_load(box_convect)
    assert_equal(str(ds), "solution-00002.pvtu")
    assert_equal(ds.dimensionality, 3)
    assert_equal(len(ds.parameters["vtu_files"]), 1)
    assert_array_equal(ds.parameters["nod_names"], _node_names)
    ic = ds.domain_center
    if can_run_ds(ds):
        dso = [None, ("sphere", (ic, (0.25, "unitary")))]
        for field in _fields:
            for dobj_name in dso:
                yield AllFieldValuesTest(ds, ("connect0", field), dobj_name)


@requires_ds(box_convect_par)
def test_box_convect_par():
    ds = data_dir_load(box_convect_par)
    assert_equal(str(ds), "solution-00001.pvtu")
    assert_equal(ds.dimensionality, 3)
    assert_equal(len(ds.parameters["vtu_files"]), 4)
    assert_array_equal(ds.parameters["nod_names"], _node_names)
    ic = ds.domain_center
    if can_run_ds(ds):
        dso = [None, ("sphere", (ic, (0.25, "unitary")))]
        for field in _fields:
            for dobj_name in dso:
                yield AllFieldValuesTest(ds, ("connect0", field), dobj_name)
