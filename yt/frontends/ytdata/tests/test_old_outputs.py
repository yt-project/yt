"""
ytdata frontend tests using enzo_tiny_cosmology



"""

import os
import shutil
import tempfile

import numpy as np

from yt.data_objects.api import create_profile
from yt.frontends.ytdata.api import (
    YTDataContainerDataset,
    YTGridDataset,
    YTNonspatialDataset,
    YTProfileDataset,
    YTSpatialPlotDataset,
)
from yt.frontends.ytdata.tests.test_outputs import (
    YTDataFieldTest,
    compare_unit_attributes,
)
from yt.testing import assert_allclose_units, assert_array_equal, requires_file, skip
from yt.units.yt_array import YTArray
from yt.utilities.answer_testing.framework import data_dir_load, requires_ds
from yt.visualization.profile_plotter import PhasePlot, ProfilePlot

enzotiny = "enzo_tiny_cosmology/DD0046/DD0046"
ytdata_dir = "ytdata_test"


@skip(reason="See https://github.com/yt-project/yt/issues/3909")
@requires_ds(enzotiny)
@requires_file(os.path.join(ytdata_dir, "DD0046_sphere.h5"))
@requires_file(os.path.join(ytdata_dir, "DD0046_cut_region.h5"))
def test_old_datacontainer_data():
    ds = data_dir_load(enzotiny)
    sphere = ds.sphere(ds.domain_center, (10, "Mpc"))
    fn = "DD0046_sphere.h5"
    full_fn = os.path.join(ytdata_dir, fn)
    sphere_ds = data_dir_load(full_fn)
    compare_unit_attributes(ds, sphere_ds)
    assert isinstance(sphere_ds, YTDataContainerDataset)
    yield YTDataFieldTest(full_fn, ("grid", "density"))
    yield YTDataFieldTest(full_fn, ("all", "particle_mass"))
    cr = ds.cut_region(sphere, ['obj[("gas", "temperature")] > 1e4'])
    fn = "DD0046_cut_region.h5"
    full_fn = os.path.join(ytdata_dir, fn)
    cr_ds = data_dir_load(full_fn)
    assert isinstance(cr_ds, YTDataContainerDataset)
    assert (cr[("gas", "temperature")] == cr_ds.data[("gas", "temperature")]).all()


@skip(reason="See https://github.com/yt-project/yt/issues/3909")
@requires_ds(enzotiny)
@requires_file(os.path.join(ytdata_dir, "DD0046_covering_grid.h5"))
@requires_file(os.path.join(ytdata_dir, "DD0046_arbitrary_grid.h5"))
@requires_file(os.path.join(ytdata_dir, "DD0046_proj_frb.h5"))
def test_old_grid_datacontainer_data():
    ds = data_dir_load(enzotiny)

    fn = "DD0046_covering_grid.h5"
    full_fn = os.path.join(ytdata_dir, fn)
    cg_ds = data_dir_load(full_fn)
    compare_unit_attributes(ds, cg_ds)
    assert isinstance(cg_ds, YTGridDataset)
    yield YTDataFieldTest(full_fn, ("grid", "density"))
    yield YTDataFieldTest(full_fn, ("all", "particle_mass"))

    fn = "DD0046_arbitrary_grid.h5"
    full_fn = os.path.join(ytdata_dir, fn)
    ag_ds = data_dir_load(full_fn)
    compare_unit_attributes(ds, ag_ds)
    assert isinstance(ag_ds, YTGridDataset)
    yield YTDataFieldTest(full_fn, ("grid", "density"))
    yield YTDataFieldTest(full_fn, ("all", "particle_mass"))

    my_proj = ds.proj("density", "x", weight_field="density")
    frb = my_proj.to_frb(1.0, (800, 800))
    fn = "DD0046_proj_frb.h5"
    full_fn = os.path.join(ytdata_dir, fn)
    frb_ds = data_dir_load(full_fn)
    assert_allclose_units(
        frb[("gas", "density")], frb_ds.data[("gas", "density")], 1e-7
    )
    compare_unit_attributes(ds, frb_ds)
    assert isinstance(frb_ds, YTGridDataset)
    yield YTDataFieldTest(full_fn, ("gas", "density"), geometric=False)


@skip(reason="See https://github.com/yt-project/yt/issues/3909")
@requires_ds(enzotiny)
@requires_file(os.path.join(ytdata_dir, "DD0046_proj.h5"))
def test_old_spatial_data():
    ds = data_dir_load(enzotiny)
    fn = "DD0046_proj.h5"
    full_fn = os.path.join(ytdata_dir, fn)
    proj_ds = data_dir_load(full_fn)
    compare_unit_attributes(ds, proj_ds)
    assert isinstance(proj_ds, YTSpatialPlotDataset)
    yield YTDataFieldTest(full_fn, ("gas", "density"), geometric=False)


@skip(reason="See https://github.com/yt-project/yt/issues/3909")
@requires_ds(enzotiny)
@requires_file(os.path.join(ytdata_dir, "DD0046_Profile1D.h5"))
@requires_file(os.path.join(ytdata_dir, "DD0046_Profile2D.h5"))
def test_old_profile_data():
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)
    ds = data_dir_load(enzotiny)
    ad = ds.all_data()
    profile_1d = create_profile(
        ad,
        ("gas", "density"),
        ("gas", "temperature"),
        weight_field=("gas", "cell_mass"),
    )
    fn = "DD0046_Profile1D.h5"
    full_fn = os.path.join(ytdata_dir, fn)
    prof_1d_ds = data_dir_load(full_fn)
    compare_unit_attributes(ds, prof_1d_ds)
    assert isinstance(prof_1d_ds, YTProfileDataset)

    for field in profile_1d.standard_deviation:
        assert_array_equal(
            profile_1d.standard_deviation[field],
            prof_1d_ds.profile.standard_deviation["data", field[1]],
        )

    p1 = ProfilePlot(
        prof_1d_ds.data,
        ("gas", "density"),
        ("gas", "temperature"),
        weight_field=("gas", "cell_mass"),
    )
    p1.save()

    yield YTDataFieldTest(full_fn, ("gas", "temperature"), geometric=False)
    yield YTDataFieldTest(full_fn, ("index", "x"), geometric=False)
    yield YTDataFieldTest(full_fn, ("gas", "density"), geometric=False)
    fn = "DD0046_Profile2D.h5"
    full_fn = os.path.join(ytdata_dir, fn)
    prof_2d_ds = data_dir_load(full_fn)
    compare_unit_attributes(ds, prof_2d_ds)
    assert isinstance(prof_2d_ds, YTProfileDataset)

    p2 = PhasePlot(
        prof_2d_ds.data,
        ("gas", "density"),
        ("gas", "temperature"),
        ("gas", "cell_mass"),
        weight_field=None,
    )
    p2.save()

    yield YTDataFieldTest(full_fn, ("gas", "density"), geometric=False)
    yield YTDataFieldTest(full_fn, ("index", "x"), geometric=False)
    yield YTDataFieldTest(full_fn, ("gas", "temperature"), geometric=False)
    yield YTDataFieldTest(full_fn, ("index", "y"), geometric=False)
    yield YTDataFieldTest(full_fn, ("gas", "cell_mass"), geometric=False)
    os.chdir(curdir)
    shutil.rmtree(tmpdir)


@skip(reason="See https://github.com/yt-project/yt/issues/3909")
@requires_ds(enzotiny)
@requires_file(os.path.join(ytdata_dir, "test_data.h5"))
@requires_file(os.path.join(ytdata_dir, "random_data.h5"))
def test_old_nonspatial_data():
    ds = data_dir_load(enzotiny)
    region = ds.box([0.25] * 3, [0.75] * 3)
    sphere = ds.sphere(ds.domain_center, (10, "Mpc"))
    my_data = {}
    my_data["region_density"] = region[("gas", "density")]
    my_data["sphere_density"] = sphere[("gas", "density")]
    fn = "test_data.h5"
    full_fn = os.path.join(ytdata_dir, fn)
    array_ds = data_dir_load(full_fn)
    compare_unit_attributes(ds, array_ds)
    assert isinstance(array_ds, YTNonspatialDataset)
    yield YTDataFieldTest(full_fn, "region_density", geometric=False)
    yield YTDataFieldTest(full_fn, "sphere_density", geometric=False)

    my_data = {"density": YTArray(np.linspace(1.0, 20.0, 10), "g/cm**3")}
    fn = "random_data.h5"
    full_fn = os.path.join(ytdata_dir, fn)
    new_ds = data_dir_load(full_fn)
    assert isinstance(new_ds, YTNonspatialDataset)
    yield YTDataFieldTest(full_fn, ("gas", "density"), geometric=False)
