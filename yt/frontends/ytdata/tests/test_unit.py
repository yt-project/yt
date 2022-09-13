import os
import shutil
import tempfile

import numpy as np

from yt.loaders import load, load_uniform_grid
from yt.testing import (
    assert_array_equal,
    assert_fname,
    fake_random_ds,
    requires_file,
    requires_module,
)
from yt.utilities.answer_testing.framework import data_dir_load
from yt.visualization.plot_window import ProjectionPlot, SlicePlot

ytdata_dir = "ytdata_test"


@requires_module("h5py")
@requires_file(os.path.join(ytdata_dir, "slice.h5"))
@requires_file(os.path.join(ytdata_dir, "proj.h5"))
@requires_file(os.path.join(ytdata_dir, "oas.h5"))
def test_old_plot_data():
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    fn = "slice.h5"
    full_fn = os.path.join(ytdata_dir, fn)
    ds_slice = data_dir_load(full_fn)
    p = SlicePlot(ds_slice, "z", ("gas", "density"))
    fn = p.save()
    assert_fname(fn[0])

    fn = "proj.h5"
    full_fn = os.path.join(ytdata_dir, fn)
    ds_proj = data_dir_load(full_fn)
    p = ProjectionPlot(ds_proj, "z", ("gas", "density"))
    fn = p.save()
    assert_fname(fn[0])

    fn = "oas.h5"
    full_fn = os.path.join(ytdata_dir, fn)
    ds_oas = data_dir_load(full_fn)
    p = SlicePlot(ds_oas, [1, 1, 1], ("gas", "density"))
    fn = p.save()
    assert_fname(fn[0])

    os.chdir(curdir)
    shutil.rmtree(tmpdir)


@requires_module("h5py")
def test_plot_data():
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)
    ds = fake_random_ds(16)

    plot = SlicePlot(ds, "z", ("gas", "density"))
    fn = plot.data_source.save_as_dataset("slice.h5")
    ds_slice = load(fn)
    p = SlicePlot(ds_slice, "z", ("gas", "density"))
    fn = p.save()
    assert_fname(fn[0])

    plot = ProjectionPlot(ds, "z", ("gas", "density"))
    fn = plot.data_source.save_as_dataset("proj.h5")
    ds_proj = load(fn)
    p = ProjectionPlot(ds_proj, "z", ("gas", "density"))
    fn = p.save()
    assert_fname(fn[0])

    plot = SlicePlot(ds, [1, 1, 1], ("gas", "density"))
    fn = plot.data_source.save_as_dataset("oas.h5")
    ds_oas = load(fn)
    p = SlicePlot(ds_oas, [1, 1, 1], ("gas", "density"))
    fn = p.save()
    assert_fname(fn[0])

    os.chdir(curdir)
    if tmpdir != ".":
        shutil.rmtree(tmpdir)


@requires_module("h5py")
def test_non_square_frb():
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    # construct an arbitrary dataset
    arr = np.arange(8.0 * 9.0 * 10.0).reshape((8, 9, 10))
    data = dict(density=(arr, "g/cm**3"))
    bbox = np.array([[-4, 4.0], [-4.5, 4.5], [-5.0, 5]])
    ds = load_uniform_grid(
        data, arr.shape, length_unit="Mpc", bbox=bbox, periodicity=(False, False, False)
    )

    # make a slice
    slc = ds.slice(axis="z", coord=ds.quan(0.0, "code_length"))
    # make a frb and save it to disk
    center = (ds.quan(0.0, "code_length"), ds.quan(0.0, "code_length"))
    xax, yax = ds.coordinates.x_axis[slc.axis], ds.coordinates.y_axis[slc.axis]
    res = [ds.domain_dimensions[xax], ds.domain_dimensions[yax]]  # = [8,9]
    width = ds.domain_right_edge[xax] - ds.domain_left_edge[xax]  # = 8 code_length
    height = ds.domain_right_edge[yax] - ds.domain_left_edge[yax]  # = 9 code_length
    frb = slc.to_frb(width=width, height=height, resolution=res, center=center)
    fname = "test_frb_roundtrip.h5"
    frb.save_as_dataset(fname, fields=[("gas", "density")])

    expected_vals = arr[:, :, 5].T
    print(
        "\nConfirmation that initial frb results are expected:",
        (expected_vals == frb[("gas", "density")].v).all(),
        "\n",
    )

    # yt-reload:
    reloaded_ds = load(fname)

    assert_array_equal(
        frb[("gas", "density")].shape, reloaded_ds.data[("gas", "density")].shape
    )
    assert_array_equal(frb[("gas", "density")], reloaded_ds.data[("gas", "density")])

    os.chdir(curdir)
    if tmpdir != ".":
        shutil.rmtree(tmpdir)
