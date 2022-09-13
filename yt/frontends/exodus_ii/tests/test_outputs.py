from yt.testing import assert_array_equal, assert_equal, requires_file
from yt.utilities.answer_testing.framework import (
    GenericArrayTest,
    data_dir_load,
    requires_ds,
)

out = "ExodusII/out.e"


@requires_file(out)
def test_out():
    ds = data_dir_load(out)
    field_list = [
        ("all", "conv_indicator"),
        ("all", "conv_marker"),
        ("all", "convected"),
        ("all", "diffused"),
        ("connect1", "conv_indicator"),
        ("connect1", "conv_marker"),
        ("connect1", "convected"),
        ("connect1", "diffused"),
        ("connect2", "conv_indicator"),
        ("connect2", "conv_marker"),
        ("connect2", "convected"),
        ("connect2", "diffused"),
    ]
    assert_equal(str(ds), "out.e")
    assert_equal(ds.dimensionality, 3)
    assert_equal(ds.current_time, 0.0)
    assert_array_equal(ds.parameters["nod_names"], ["convected", "diffused"])
    assert_equal(ds.parameters["num_meshes"], 2)
    assert_array_equal(ds.field_list, field_list)


out_s002 = "ExodusII/out.e-s002"


@requires_file(out_s002)
def test_out002():
    ds = data_dir_load(out_s002)
    field_list = [
        ("all", "conv_indicator"),
        ("all", "conv_marker"),
        ("all", "convected"),
        ("all", "diffused"),
        ("connect1", "conv_indicator"),
        ("connect1", "conv_marker"),
        ("connect1", "convected"),
        ("connect1", "diffused"),
        ("connect2", "conv_indicator"),
        ("connect2", "conv_marker"),
        ("connect2", "convected"),
        ("connect2", "diffused"),
    ]
    assert_equal(str(ds), "out.e-s002")
    assert_equal(ds.dimensionality, 3)
    assert_equal(ds.current_time, 2.0)
    assert_array_equal(ds.field_list, field_list)


gold = "ExodusII/gold.e"


@requires_file(gold)
def test_gold():
    ds = data_dir_load(gold)
    field_list = [("all", "forced"), ("connect1", "forced")]
    assert_equal(str(ds), "gold.e")
    assert_array_equal(ds.field_list, field_list)


big_data = "MOOSE_sample_data/mps_out.e"


@requires_ds(big_data)
def test_displacement_fields():
    displacement_dicts = [
        {"connect2": (5.0, [0.0, 0.0, 0.0])},
        {"connect1": (1.0, [1.0, 2.0, 3.0]), "connect2": (0.0, [0.0, 0.0, 0.0])},
    ]
    for disp in displacement_dicts:
        ds = data_dir_load(big_data, displacements=disp)
        for mesh in ds.index.meshes:

            def array_func():
                return mesh.connectivity_coords  # noqa: B023

            yield GenericArrayTest(ds, array_func, 12)
