from yt.testing import assert_array_equal, assert_raises, fake_amr_ds, fake_random_ds
from yt.units import cm


def test_box_creation():

    # test that creating a region with left and right edge
    # with units works
    ds = fake_random_ds(32, length_unit=2)
    reg = ds.box([0, 0, 0] * cm, [2, 2, 2] * cm)
    dens_units = reg["density"]

    reg = ds.box([0, 0, 0], [1, 1, 1])
    dens_no_units = reg["density"]

    assert_array_equal(dens_units, dens_no_units)


def test_max_level_min_level_semantics():
    ds = fake_amr_ds()
    ad = ds.all_data()
    assert ad["grid_level"].max() == 4
    ad.max_level = 2
    assert ad["grid_level"].max() == 2
    ad.max_level = 8
    assert ad["grid_level"].max() == 4
    ad.min_level = 2
    assert ad["grid_level"].min() == 2
    ad.min_level = 0
    assert ad["grid_level"].min() == 0


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
