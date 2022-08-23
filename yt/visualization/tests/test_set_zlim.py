import numpy as np
import numpy.testing as npt
import pytest

from yt._maintenance.deprecation import VisibleDeprecationWarning
from yt.testing import fake_amr_ds
from yt.visualization.api import SlicePlot


def test_float_vmin_then_set_unit():
    # this test doesn't represent how users should interact with plot containers
    # in particular it uses the `_setup_plots()` private method, as a quick way to
    # create a plot without having to make it an answer test
    field = ("gas", "density")
    ds = fake_amr_ds(fields=[field], units=["g/cm**3"])

    p = SlicePlot(ds, "x", field)
    p.set_buff_size(16)

    p._setup_plots()
    cb = p.plots[field].image.colorbar
    raw_lims = np.array((cb.vmin, cb.vmax))
    desired_lims = raw_lims.copy()
    desired_lims[0] = 1e-2

    p.set_zlim(field, zmin=desired_lims[0])

    p._setup_plots()
    cb = p.plots[field].image.colorbar
    new_lims = np.array((cb.vmin, cb.vmax))
    npt.assert_almost_equal(new_lims, desired_lims)

    # 1 g/cm**3 == 1000 kg/m**3
    p.set_unit(field, "kg/m**3")
    p._setup_plots()

    cb = p.plots[field].image.colorbar
    new_lims = np.array((cb.vmin, cb.vmax))
    npt.assert_almost_equal(new_lims, 1000 * desired_lims)


def test_set_unit_then_float_vmin():
    field = ("gas", "density")
    ds = fake_amr_ds(fields=[field], units=["g/cm**3"])

    p = SlicePlot(ds, "x", field)
    p.set_buff_size(16)

    p.set_unit(field, "kg/m**3")
    p.set_zlim(field, zmin=1)
    p._setup_plots()
    cb = p.plots[field].image.colorbar
    assert cb.vmin == 1.0


def test_reset_zlim():
    field = ("gas", "density")
    ds = fake_amr_ds(fields=[field], units=["g/cm**3"])

    p = SlicePlot(ds, "x", field)
    p.set_buff_size(16)

    p._setup_plots()
    cb = p.plots[field].image.colorbar
    raw_lims = np.array((cb.vmin, cb.vmax))

    # set a new zmin value
    delta = np.diff(raw_lims)[0]
    p.set_zlim(field, zmin=raw_lims[0] + delta / 2)

    # passing "min" should restore default limit
    p.set_zlim(field, zmin="min")
    p._setup_plots()

    cb = p.plots[field].image.colorbar
    new_lims = np.array((cb.vmin, cb.vmax))
    npt.assert_array_equal(new_lims, raw_lims)


def test_set_dynamic_range_with_vmin():
    field = ("gas", "density")
    ds = fake_amr_ds(fields=[field], units=["g/cm**3"])

    p = SlicePlot(ds, "x", field)
    p.set_buff_size(16)

    zmin = 1e-2
    p.set_zlim(field, zmin=zmin, dynamic_range=2)

    p._setup_plots()
    cb = p.plots[field].image.colorbar
    new_lims = np.array((cb.vmin, cb.vmax))
    npt.assert_almost_equal(new_lims, (zmin, 2 * zmin))


def test_set_dynamic_range_with_vmax():
    field = ("gas", "density")
    ds = fake_amr_ds(fields=[field], units=["g/cm**3"])

    p = SlicePlot(ds, "x", field)
    p.set_buff_size(16)

    zmax = 1
    p.set_zlim(field, zmax=zmax, dynamic_range=2)

    p._setup_plots()
    cb = p.plots[field].image.colorbar
    new_lims = np.array((cb.vmin, cb.vmax))
    npt.assert_almost_equal(new_lims, (zmax / 2, zmax))


def test_set_dynamic_range_with_min():
    field = ("gas", "density")
    ds = fake_amr_ds(fields=[field], units=["g/cm**3"])

    p = SlicePlot(ds, "x", field)
    p.set_buff_size(16)

    p._setup_plots()
    cb = p.plots[field].image.colorbar
    vmin = cb.vmin

    p.set_zlim(field, zmin="min", dynamic_range=2)

    p._setup_plots()
    cb = p.plots[field].image.colorbar
    new_lims = np.array((cb.vmin, cb.vmax))
    npt.assert_almost_equal(new_lims, (vmin, 2 * vmin))


def test_set_dynamic_range_with_None():
    field = ("gas", "density")
    ds = fake_amr_ds(fields=[field], units=["g/cm**3"])

    p = SlicePlot(ds, "x", field)
    p.set_buff_size(16)

    p._setup_plots()
    cb = p.plots[field].image.colorbar
    vmin = cb.vmin

    with pytest.raises(
        VisibleDeprecationWarning, match="Passing `zmin=None` explicitly is deprecated"
    ):
        p.set_zlim(field, zmin=None, dynamic_range=2)

        p._setup_plots()
        cb = p.plots[field].image.colorbar
        new_lims = np.array((cb.vmin, cb.vmax))
        npt.assert_almost_equal(new_lims, (vmin, 2 * vmin))
