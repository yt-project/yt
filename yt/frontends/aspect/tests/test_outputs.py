from yt.testing import assert_array_equal, assert_equal, requires_file
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


def common_assertions(ds, solution_step, nproc):
    assert_equal(str(ds), f"solution-0000{solution_step}.pvtu")
    assert_equal(ds.dimensionality, 3)
    assert_equal(len(ds.parameters["vtu_files"]), nproc)
    assert_array_equal(ds.parameters["nod_names"], _node_names)


@requires_file(box_convect)
def test_box_convect():
    ds = data_dir_load(box_convect)
    common_assertions(ds, 2, 1)


@requires_file(box_convect_par)
def test_box_convect_par():
    ds = data_dir_load(box_convect_par)
    common_assertions(ds, 1, 4)


@requires_ds(box_convect)
def test_box_convect_values():
    ds = data_dir_load(box_convect)
    ic = ds.domain_center
    if can_run_ds(ds):
        dso = [None, ("sphere", (ic, (0.25, "unitary")))]
        for field in _fields:
            for dobj_name in dso:
                yield AllFieldValuesTest(ds, field, dobj_name)


@requires_ds(box_convect_par)
def test_box_convect_par_values():
    ds = data_dir_load(box_convect_par)
    ic = ds.domain_center
    if can_run_ds(ds):
        dso = [None, ("sphere", (ic, (0.25, "unitary")))]
        for field in _fields:
            for dobj_name in dso:
                yield AllFieldValuesTest(ds, field, dobj_name)
