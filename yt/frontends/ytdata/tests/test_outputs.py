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
    save_as_dataset,
)
from yt.loaders import load
from yt.testing import assert_allclose_units, assert_array_equal, assert_equal
from yt.units.yt_array import YTArray, YTQuantity
from yt.utilities.answer_testing.framework import (
    AnswerTestingTest,
    data_dir_load,
    requires_ds,
)
from yt.visualization.profile_plotter import PhasePlot, ProfilePlot


def make_tempdir():
    if int(os.environ.get("GENERATE_YTDATA", 0)):
        return "."
    else:
        return tempfile.mkdtemp()


def compare_unit_attributes(ds1, ds2):
    r"""
    Checks to make sure that the length, mass, time, velocity, and
    magnetic units are the same for two different dataset objects.
    """
    attrs = ("length_unit", "mass_unit", "time_unit", "velocity_unit", "magnetic_unit")
    for attr in attrs:
        u1 = getattr(ds1, attr, None)
        u2 = getattr(ds2, attr, None)
        assert u1 == u2


class YTDataFieldTest(AnswerTestingTest):
    _type_name = "YTDataTest"
    _attrs = ("field_name",)

    def __init__(self, ds_fn, field, decimals=10, geometric=True):
        super().__init__(ds_fn)
        self.field = field
        if isinstance(field, tuple) and len(field) == 2:
            self.field_name = field[1]
        else:
            self.field_name = field
        self.decimals = decimals
        self.geometric = geometric

    def run(self):
        if self.geometric:
            obj = self.ds.all_data()
        else:
            obj = self.ds.data
        num_e = obj[self.field].size
        avg = obj[self.field].mean()
        return np.array([num_e, avg])

    def compare(self, new_result, old_result):
        err_msg = f"YTData field values for {self.field} not equal."
        if self.decimals is None:
            assert_equal(new_result, old_result, err_msg=err_msg, verbose=True)
        else:
            assert_allclose_units(
                new_result,
                old_result,
                10.0 ** (-self.decimals),
                err_msg=err_msg,
                verbose=True,
            )


enzotiny = "enzo_tiny_cosmology/DD0046/DD0046"


@requires_ds(enzotiny)
def test_datacontainer_data():
    tmpdir = make_tempdir()
    curdir = os.getcwd()
    os.chdir(tmpdir)
    ds = data_dir_load(enzotiny)
    sphere = ds.sphere(ds.domain_center, (10, "Mpc"))
    fn = sphere.save_as_dataset(fields=[("gas", "density"), ("all", "particle_mass")])
    full_fn = os.path.join(tmpdir, fn)
    sphere_ds = load(full_fn)
    compare_unit_attributes(ds, sphere_ds)
    assert isinstance(sphere_ds, YTDataContainerDataset)
    yield YTDataFieldTest(full_fn, ("grid", "density"))
    yield YTDataFieldTest(full_fn, ("all", "particle_mass"))
    cr = ds.cut_region(sphere, ['obj[("gas", "temperature")] > 1e4'])
    fn = cr.save_as_dataset(fields=[("gas", "temperature")])
    full_fn = os.path.join(tmpdir, fn)
    cr_ds = load(full_fn)
    assert isinstance(cr_ds, YTDataContainerDataset)
    assert (cr[("gas", "temperature")] == cr_ds.data[("gas", "temperature")]).all()
    os.chdir(curdir)
    if tmpdir != ".":
        shutil.rmtree(tmpdir)


@requires_ds(enzotiny)
def test_grid_datacontainer_data():
    tmpdir = make_tempdir()
    curdir = os.getcwd()
    os.chdir(tmpdir)
    ds = data_dir_load(enzotiny)

    cg = ds.covering_grid(level=0, left_edge=[0.25] * 3, dims=[16] * 3)
    fn = cg.save_as_dataset(
        fields=[
            ("gas", "density"),
            ("all", "particle_mass"),
            ("all", "particle_position"),
        ]
    )
    full_fn = os.path.join(tmpdir, fn)
    cg_ds = load(full_fn)
    compare_unit_attributes(ds, cg_ds)
    assert isinstance(cg_ds, YTGridDataset)
    assert (
        cg["all", "particle_position"].shape
        == cg_ds.r["all", "particle_position"].shape
    )
    yield YTDataFieldTest(full_fn, ("grid", "density"))
    yield YTDataFieldTest(full_fn, ("all", "particle_mass"))

    ag = ds.arbitrary_grid(left_edge=[0.25] * 3, right_edge=[0.75] * 3, dims=[16] * 3)
    fn = ag.save_as_dataset(fields=[("gas", "density"), ("all", "particle_mass")])
    full_fn = os.path.join(tmpdir, fn)
    ag_ds = load(full_fn)
    compare_unit_attributes(ds, ag_ds)
    assert isinstance(ag_ds, YTGridDataset)
    yield YTDataFieldTest(full_fn, ("grid", "density"))
    yield YTDataFieldTest(full_fn, ("all", "particle_mass"))

    my_proj = ds.proj(("gas", "density"), "x", weight_field=("gas", "density"))
    frb = my_proj.to_frb(1.0, (800, 800))
    fn = frb.save_as_dataset(fields=[("gas", "density")])
    frb_ds = load(fn)
    assert_array_equal(frb[("gas", "density")], frb_ds.data[("gas", "density")])
    compare_unit_attributes(ds, frb_ds)
    assert isinstance(frb_ds, YTGridDataset)
    yield YTDataFieldTest(full_fn, ("grid", "density"), geometric=False)
    os.chdir(curdir)
    if tmpdir != ".":
        shutil.rmtree(tmpdir)


@requires_ds(enzotiny)
def test_spatial_data():
    tmpdir = make_tempdir()
    curdir = os.getcwd()
    os.chdir(tmpdir)
    ds = data_dir_load(enzotiny)
    proj = ds.proj(("gas", "density"), "x", weight_field=("gas", "density"))
    fn = proj.save_as_dataset()
    full_fn = os.path.join(tmpdir, fn)
    proj_ds = load(full_fn)
    compare_unit_attributes(ds, proj_ds)
    assert isinstance(proj_ds, YTSpatialPlotDataset)
    yield YTDataFieldTest(full_fn, ("grid", "density"), geometric=False)
    os.chdir(curdir)
    if tmpdir != ".":
        shutil.rmtree(tmpdir)


@requires_ds(enzotiny)
def test_profile_data():
    tmpdir = make_tempdir()
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
    fn = profile_1d.save_as_dataset()
    full_fn = os.path.join(tmpdir, fn)
    prof_1d_ds = load(full_fn)
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

    yield YTDataFieldTest(full_fn, ("data", "temperature"), geometric=False)
    yield YTDataFieldTest(full_fn, ("data", "x"), geometric=False)
    yield YTDataFieldTest(full_fn, ("data", "density"), geometric=False)
    profile_2d = create_profile(
        ad,
        [("gas", "density"), ("gas", "temperature")],
        ("gas", "cell_mass"),
        weight_field=None,
        n_bins=(128, 128),
    )
    fn = profile_2d.save_as_dataset()
    full_fn = os.path.join(tmpdir, fn)
    prof_2d_ds = load(full_fn)
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

    yield YTDataFieldTest(full_fn, ("data", "density"), geometric=False)
    yield YTDataFieldTest(full_fn, ("data", "x"), geometric=False)
    yield YTDataFieldTest(full_fn, ("data", "temperature"), geometric=False)
    yield YTDataFieldTest(full_fn, ("data", "y"), geometric=False)
    yield YTDataFieldTest(full_fn, ("data", "cell_mass"), geometric=False)
    os.chdir(curdir)
    if tmpdir != ".":
        shutil.rmtree(tmpdir)


@requires_ds(enzotiny)
def test_nonspatial_data():
    tmpdir = make_tempdir()
    curdir = os.getcwd()
    os.chdir(tmpdir)
    ds = data_dir_load(enzotiny)
    region = ds.box([0.25] * 3, [0.75] * 3)
    sphere = ds.sphere(ds.domain_center, (10, "Mpc"))
    my_data = {}
    my_data["region_density"] = region[("gas", "density")]
    my_data["sphere_density"] = sphere[("gas", "density")]
    fn = "test_data.h5"
    save_as_dataset(ds, fn, my_data)
    full_fn = os.path.join(tmpdir, fn)
    array_ds = load(full_fn)
    compare_unit_attributes(ds, array_ds)
    assert isinstance(array_ds, YTNonspatialDataset)
    yield YTDataFieldTest(full_fn, "region_density", geometric=False)
    yield YTDataFieldTest(full_fn, "sphere_density", geometric=False)

    my_data = {"density": YTArray(np.linspace(1.0, 20.0, 10), "g/cm**3")}
    fake_ds = {"current_time": YTQuantity(10, "Myr")}
    fn = "random_data.h5"
    save_as_dataset(fake_ds, fn, my_data)
    full_fn = os.path.join(tmpdir, fn)
    new_ds = load(full_fn)
    assert isinstance(new_ds, YTNonspatialDataset)
    yield YTDataFieldTest(full_fn, ("data", "density"), geometric=False)
    os.chdir(curdir)
    if tmpdir != ".":
        shutil.rmtree(tmpdir)
