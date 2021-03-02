from yt.testing import assert_array_equal, assert_equal, requires_file
from yt.utilities.answer_testing.framework import data_dir_load

_node_names = ("velocity", "p", "T")
_fields = ("T", "p", "velocity_x", "velocity_z", "velocity_y")
box_convect = "aspect/output_convection_box_3d/nproc_1/solution/solution-00002.pvtu"
box_convect_par = "aspect/output_convection_box_3d/nproc_4/solution/solution-00002.pvtu"


def common_assertions(ds, nproc):
    assert_equal(str(ds), "solution-00002.pvtu")
    assert_equal(ds.dimensionality, 3)
    assert_equal(len(ds.parameters["vtu_files"]), nproc)
    assert_array_equal(ds.parameters["nod_names"], _node_names)


@requires_file(box_convect)
def test_box_convect():
    ds = data_dir_load(box_convect)
    common_assertions(ds, 1)


@requires_file(box_convect_par)
def test_box_convect_par():
    ds = data_dir_load(box_convect_par)
    common_assertions(ds, 4)
