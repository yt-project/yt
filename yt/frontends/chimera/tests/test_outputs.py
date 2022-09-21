"""
Chimera frontend tests

"""

import numpy as np

from yt.testing import (
    assert_almost_equal,
    assert_array_equal,
    assert_equal,
    requires_file,
)
from yt.utilities.answer_testing.framework import (
    GenericArrayTest,
    data_dir_load,
    requires_ds,
)

Two_D = "F37_80/chimera_00001_grid_1_01.h5"


@requires_ds(Two_D)
def test_2D():
    ds = data_dir_load(Two_D)
    _fields = [
        ("chimera", "a_nuc_rep_c"),
        ("chimera", "abar"),
        ("chimera", "ar36"),
        ("chimera", "be_nuc_rep_c"),
        ("chimera", "c12"),
        ("chimera", "ca40"),
        ("chimera", "cr48"),
        ("chimera", "dudt_nu"),
        ("chimera", "dudt_nuc"),
        ("chimera", "e_book"),
        ("chimera", "e_int"),
        ("chimera", "e_rms_1"),
        ("chimera", "e_rms_2"),
        ("chimera", "e_rms_3"),
        ("chimera", "e_rms_4"),
        ("chimera", "entropy"),
        ("chimera", "fe52"),
        ("chimera", "fe56"),
        ("chimera", "grav_x_c"),
        ("chimera", "grav_y_c"),
        ("chimera", "grav_z_c"),
        ("chimera", "he4"),
        ("chimera", "lumin_1"),
        ("chimera", "lumin_2"),
        ("chimera", "lumin_3"),
        ("chimera", "lumin_4"),
        ("chimera", "mg24"),
        ("chimera", "n"),
        ("chimera", "ne20"),
        ("chimera", "ni56"),
        ("chimera", "nse_c"),
        ("chimera", "num_lumin_1"),
        ("chimera", "num_lumin_2"),
        ("chimera", "num_lumin_3"),
        ("chimera", "num_lumin_4"),
        ("chimera", "o16"),
        ("chimera", "p"),
        ("chimera", "press"),
        ("chimera", "rho_c"),
        ("chimera", "s32"),
        ("chimera", "si28"),
        ("chimera", "t_c"),
        ("chimera", "ti44"),
        ("chimera", "u_c"),
        ("chimera", "v_c"),
        ("chimera", "v_csound"),
        ("chimera", "wBVMD"),
        ("chimera", "w_c"),
        ("chimera", "ye_c"),
        ("chimera", "ylep"),
        ("chimera", "z_nuc_rep_c"),
        ("chimera", "zn60"),
    ]
    assert_equal(str(ds), "chimera_00001_grid_1_01.h5")
    assert_equal(ds.geometry, "spherical")  # Geometry
    assert_almost_equal(
        ds.domain_right_edge,
        ds.arr([1.0116509e10 + 100, 3.14159265e00, 6.28318531e00], "code_length"),
    )  # domain edge
    assert_array_equal(
        ds.domain_left_edge, ds.arr([0.0, 0.0, 0.0], "code_length")
    )  # domain edge
    assert_array_equal(ds.domain_dimensions, np.array([722, 240, 1]))  # Dimensions
    assert_array_equal(ds.field_list, _fields)

    def field_func(field):
        min = dd[field].min()
        max = dd[field].max()
        avg = np.mean(dd[field])
        size = dd[field].size
        return [min, max, avg, size]

    dd = ds.all_data()
    for field in _fields:

        if field != ("chimera", "shock"):
            yield GenericArrayTest(ds, field_func, args=[field])


Three_D = "C15-3D-3deg/chimera_002715000_grid_1_01.h5"


@requires_ds(Three_D)
def test_3D():
    ds = data_dir_load(Three_D)
    _fields = [
        ("chimera", "a_nuc_rep_c"),
        ("chimera", "abar"),
        ("chimera", "ar36"),
        ("chimera", "be_nuc_rep_c"),
        ("chimera", "c12"),
        ("chimera", "ca40"),
        ("chimera", "cr48"),
        ("chimera", "dudt_nu"),
        ("chimera", "dudt_nuc"),
        ("chimera", "e_book"),
        ("chimera", "e_int"),
        ("chimera", "e_rms_1"),
        ("chimera", "e_rms_2"),
        ("chimera", "e_rms_3"),
        ("chimera", "e_rms_4"),
        ("chimera", "entropy"),
        ("chimera", "fe52"),
        ("chimera", "grav_x_c"),
        ("chimera", "grav_y_c"),
        ("chimera", "grav_z_c"),
        ("chimera", "he4"),
        ("chimera", "lumin_1"),
        ("chimera", "lumin_2"),
        ("chimera", "lumin_3"),
        ("chimera", "lumin_4"),
        ("chimera", "mg24"),
        ("chimera", "n"),
        ("chimera", "ne20"),
        ("chimera", "ni56"),
        ("chimera", "nse_c"),
        ("chimera", "num_lumin_1"),
        ("chimera", "num_lumin_2"),
        ("chimera", "num_lumin_3"),
        ("chimera", "num_lumin_4"),
        ("chimera", "o16"),
        ("chimera", "p"),
        ("chimera", "press"),
        ("chimera", "rho_c"),
        ("chimera", "s32"),
        ("chimera", "si28"),
        ("chimera", "t_c"),
        ("chimera", "ti44"),
        ("chimera", "u_c"),
        ("chimera", "v_c"),
        ("chimera", "v_csound"),
        ("chimera", "wBVMD"),
        ("chimera", "w_c"),
        ("chimera", "ye_c"),
        ("chimera", "ylep"),
        ("chimera", "z_nuc_rep_c"),
        ("chimera", "zn60"),
    ]
    assert_equal(str(ds), "chimera_002715000_grid_1_01.h5")
    assert_equal(ds.geometry, "spherical")  # Geometry
    assert_almost_equal(
        ds.domain_right_edge,
        ds.arr(
            [1.06500257e09 - 1.03818333, 3.14159265e00, 6.2831853e00], "code_length"
        ),
    )  # Domain edge
    assert_array_equal(ds.domain_left_edge, [0.0, 0.0, 0.0])  # Domain edge
    assert_array_equal(ds.domain_dimensions, [542, 60, 135])  # Dimensions
    assert_array_equal(ds.field_list, _fields)

    def field_func(field):
        min = dd[field].min()
        max = dd[field].max()
        avg = np.mean(dd[field])
        size = dd[field].size
        return [min, max, avg, size]

    dd = ds.all_data()
    for field in _fields:

        if field != ("chimera", "shock"):
            yield GenericArrayTest(ds, field_func, args=[field])


@requires_file(Three_D)
def test_multimesh():  # Tests that the multimesh system for 3D data has been created correctly
    ds = data_dir_load(Three_D)
    assert_equal(len(ds.index.meshes), 45)
    for i in range(44):
        assert_almost_equal(
            ds.index.meshes[i + 1].connectivity_coords
            - ds.index.meshes[i].connectivity_coords,
            np.tile([0.0, 0.0, 0.13962634015954636], (132004, 1)),
        )  # Tests that each mesh is an identically shaped wedge, incrememnted in Phi.
        assert_array_equal(
            ds.index.meshes[i + 1].connectivity_indices,
            ds.index.meshes[i].connectivity_indices,
        )  # Checks Connectivity array is identical for all meshes.
