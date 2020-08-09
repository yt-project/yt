import os
import shutil
import tempfile

from yt.loaders import load
from yt.testing import assert_fname, fake_random_ds, requires_file, requires_module
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
    p = SlicePlot(ds_slice, "z", "density")
    fn = p.save()
    assert_fname(fn[0])

    fn = "proj.h5"
    full_fn = os.path.join(ytdata_dir, fn)
    ds_proj = data_dir_load(full_fn)
    p = ProjectionPlot(ds_proj, "z", "density")
    fn = p.save()
    assert_fname(fn[0])

    fn = "oas.h5"
    full_fn = os.path.join(ytdata_dir, fn)
    ds_oas = data_dir_load(full_fn)
    p = SlicePlot(ds_oas, [1, 1, 1], "density")
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

    plot = SlicePlot(ds, "z", "density")
    fn = plot.data_source.save_as_dataset("slice.h5")
    ds_slice = load(fn)
    p = SlicePlot(ds_slice, "z", "density")
    fn = p.save()
    assert_fname(fn[0])

    plot = ProjectionPlot(ds, "z", "density")
    fn = plot.data_source.save_as_dataset("proj.h5")
    ds_proj = load(fn)
    p = ProjectionPlot(ds_proj, "z", "density")
    fn = p.save()
    assert_fname(fn[0])

    plot = SlicePlot(ds, [1, 1, 1], "density")
    fn = plot.data_source.save_as_dataset("oas.h5")
    ds_oas = load(fn)
    p = SlicePlot(ds_oas, [1, 1, 1], "density")
    fn = p.save()
    assert_fname(fn[0])

    os.chdir(curdir)
    if tmpdir != ".":
        shutil.rmtree(tmpdir)
