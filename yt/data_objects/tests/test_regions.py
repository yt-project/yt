import numpy as np
from numpy.testing import assert_array_equal, assert_equal, assert_raises

import yt
from yt.testing import fake_amr_ds, fake_random_ds
from yt.units import cm


def test_box_creation():
    # test that creating a region with left and right edge
    # with units works
    ds = fake_random_ds(32, length_unit=2)
    reg = ds.box([0, 0, 0] * cm, [2, 2, 2] * cm)
    dens_units = reg["gas", "density"]

    reg = ds.box([0, 0, 0], [1, 1, 1])
    dens_no_units = reg["gas", "density"]

    assert_array_equal(dens_units, dens_no_units)


def test_box_physical_unit_edges_match_code_length():
    # https://github.com/yt-project/yt/issues/5496
    # A subset box in physical units must select the same particles as the
    # equivalent code_length box. The existing test_box_creation uses the
    # full domain, which hides the mismatch.
    positions = np.array([0.10, 0.45, 0.50, 0.55, 0.90])
    data = {
        "particle_position_x": positions,
        "particle_position_y": positions,
        "particle_position_z": positions,
        "particle_mass": np.ones(positions.size),
    }
    ds = yt.load_particles(data, length_unit=(10.0, "Mpc"))
    center_code = ds.arr([0.5, 0.5, 0.5], "code_length")
    center_mpc = center_code.to("Mpc")
    half_width_code = ds.quan(0.1, "code_length")
    half_width_mpc = half_width_code.to("Mpc")

    box_code = ds.box(center_code - half_width_code, center_code + half_width_code)
    box_mpc = ds.box(center_mpc - half_width_mpc, center_mpc + half_width_mpc)
    n_code = box_code["all", "particle_mass"].size
    n_mpc = box_mpc["all", "particle_mass"].size
    assert_equal(n_code, 3)
    assert_equal(n_mpc, n_code)

    reg_code = ds.region(center_code, center_code - half_width_code, center_code + half_width_code)
    reg_mpc = ds.region(center_mpc, center_mpc - half_width_mpc, center_mpc + half_width_mpc)
    assert_equal(reg_mpc["all", "particle_mass"].size, reg_code["all", "particle_mass"].size)

    sph_code = ds.sphere(center_code, half_width_code)
    sph_mpc = ds.sphere(center_mpc, half_width_mpc)
    assert_equal(sph_code["all", "particle_mass"].size, 3)
    assert_equal(sph_mpc["all", "particle_mass"].size, 3)


def test_max_level_min_level_semantics():
    ds = fake_amr_ds()
    ad = ds.all_data()
    assert ad["index", "grid_level"].max() == 4
    ad.max_level = 2
    assert ad["index", "grid_level"].max() == 2
    ad.max_level = 8
    assert ad["index", "grid_level"].max() == 4
    ad.min_level = 2
    assert ad["index", "grid_level"].min() == 2
    ad.min_level = 0
    assert ad["index", "grid_level"].min() == 0


def test_ellipsis_selection():
    ds = fake_amr_ds()
    reg = ds.r[:, :, :]
    ereg = ds.r[...]
    assert_array_equal(reg.fwidth, ereg.fwidth)

    reg = ds.r[(0.5, "cm"), :, :]
    ereg = ds.r[(0.5, "cm"), ...]
    assert_array_equal(reg.fwidth, ereg.fwidth)

    reg = ds.r[:, :, (0.5, "cm")]
    ereg = ds.r[..., (0.5, "cm")]
    assert_array_equal(reg.fwidth, ereg.fwidth)

    reg = ds.r[:, :, (0.5, "cm")]
    ereg = ds.r[..., (0.5, "cm")]
    assert_array_equal(reg.fwidth, ereg.fwidth)

    assert_raises(IndexError, ds.r.__getitem__, (..., (0.5, "cm"), ...))


# this test will fail until "arbitrary_grid" selector is implemented for 2D datasets
# see issue https://github.com/yt-project/yt/issues/3437
"""
from yt.utilities.answer_testing.framework import data_dir_load, requires_ds

@requires_ds("castro_sedov_2d_cyl_in_cart_plt00150")
def test_complex_slicing_2D_consistency():

    # see https://github.com/yt-project/yt/issues/3429
    ds = data_dir_load("castro_sedov_2d_cyl_in_cart_plt00150")

    reg = ds.r[0.1:0.2:8j, :]
    reg["gas", "density"]
    reg = ds.r[:, 1:2:8j]
    reg["gas", "density"]
    reg = ds.r[0.1:0.2:8j, 1:2:8j]
    reg["gas", "density"]
"""
