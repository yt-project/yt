"""
ytdata frontend tests using enzo_tiny_cosmology



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.convenience import \
    load
from yt.frontends.ytdata.api import \
    YTDataContainerDataset, \
    YTSpatialPlotDataset, \
    YTGridDataset, \
    YTNonspatialDataset, \
    YTProfileDataset, \
    save_as_dataset
from yt.testing import \
    assert_allclose_units, \
    assert_equal
from yt.utilities.answer_testing.framework import \
    requires_ds, \
    data_dir_load, \
    AnswerTestingTest
from yt.units.yt_array import \
    YTArray, \
    YTQuantity
from yt.data_objects.api import \
    create_profile
import numpy as np
import tempfile
import os
import shutil

class YTDataFieldTest(AnswerTestingTest):
    _type_name = "YTDataTest"
    _attrs = ("field_name", )

    def __init__(self, ds_fn, field, decimals = 10,
                 geometric=True):
        super(YTDataFieldTest, self).__init__(ds_fn)
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
        num_e = obj[field].size
        avg = obj[field].mean()
        return np.array([num_e, avg])

    def compare(self, new_result, old_result):
        err_msg = "YTData field values for %s not equal." % \
          (self.field,)
        if self.decimals is None:
            assert_equal(new_result, old_result,
                         err_msg=err_msg, verbose=True)
        else:
            assert_allclose_units(new_result, old_result, 
                                  10.**(-self.decimals),
                                  err_msg=err_msg, verbose=True)

enzotiny = "enzo_tiny_cosmology/DD0046/DD0046"
@requires_ds(enzotiny)
def test_datacontainer_data():
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)
    ds = data_dir_load(enzotiny)
    sphere = ds.sphere(ds.domain_center, (10, "Mpc"))
    fn = sphere.save_as_dataset(fields=["density", "particle_mass"])
    sphere_ds = load(fn)
    assert isinstance(sphere_ds, YTDataContainerDataset)
    yield YTDataFieldTest(fn, ("grid", "density"))
    yield YTDataFieldTest(fn, ("all", "particle_mass"))
    os.chdir(curdir)
    shutil.rmtree(tmpdir)

@requires_ds(enzotiny)
def test_grid_datacontainer_data():
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)
    ds = data_dir_load(enzotiny)
    cg = ds.covering_grid(level=0, left_edge=[0.25]*3, dims=[16]*3)
    fn = cg.save_as_dataset(fields=["density", "particle_mass"])
    cg_ds = load(fn)
    assert isinstance(cg_ds, YTGridDataset)
    yield YTDataFieldTest(fn, ("grid", "density"))
    yield YTDataFieldTest(fn, ("all", "particle_mass"))

    my_proj = ds.proj("density", "x", weight_field="density")
    frb = my_proj.to_frb(1.0, (800, 800))
    fn = frb.save_as_dataset(fields=["density"])
    frb_ds = load(fn)
    assert isinstance(frb_ds, YTGridDataset)
    yield YTDataFieldTest(fn, "density", geometric=False)
    os.chdir(curdir)
    shutil.rmtree(tmpdir)

@requires_ds(enzotiny)
def test_spatial_data():
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)
    ds = data_dir_load(enzotiny)
    proj = ds.proj("density", "x", weight_field="density")
    fn = proj.save_as_dataset()
    proj_ds = load(fn)
    assert isinstance(proj_ds, YTSpatialPlotDataset)
    yield YTDataFieldTest(fn, ("grid", "density"), geometric=False)
    os.chdir(curdir)
    shutil.rmtree(tmpdir)

@requires_ds(enzotiny)
def test_profile_data():
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)
    ds = data_dir_load(enzotiny)
    ad = ds.all_data()
    profile_1d = create_profile(ad, "density", "temperature",
                                weight_field="cell_mass")
    fn = profile_1d.save_as_dataset()
    prof_1d_ds = load(fn)
    assert isinstance(prof_1d_ds, YTProfileDataset)
    yield YTDataFieldTest(fn, "temperature", geometric=False)
    yield YTDataFieldTest(fn, "x", geometric=False)
    yield YTDataFieldTest(fn, "density", geometric=False)

    profile_2d = create_profile(ad, ["density", "temperature"],
                               "cell_mass", weight_field=None,
                               n_bins=(128, 128))
    fn = profile_2d.save_as_dataset()
    prof_2d_ds = load(fn)
    assert isinstance(prof_2d_ds, YTProfileDataset)
    yield YTDataFieldTest(fn, "density", geometric=False)
    yield YTDataFieldTest(fn, "x", geometric=False)
    yield YTDataFieldTest(fn, "temperature", geometric=False)
    yield YTDataFieldTest(fn, "y", geometric=False)
    yield YTDataFieldTest(fn, "cell_mass", geometric=False)
    os.chdir(curdir)
    shutil.rmtree(tmpdir)

@requires_ds(enzotiny)
def test_nonspatial_data():
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)
    ds = data_dir_load(enzotiny)
    region = ds.box([0.25]*3, [0.75]*3)
    sphere = ds.sphere(ds.domain_center, (10, "Mpc"))
    my_data = {}
    my_data["region_density"] = region["density"]
    my_data["sphere_density"] = sphere["density"]
    fn = save_as_dataset(ds, "test_data.h5", my_data)
    array_ds = load(fn)
    assert isinstance(array_ds, YTNonspatialDataset)
    yield YTDataFieldTest(fn, "region_density", geometric=False)
    yield YTDataFieldTest(fn, "sphere_density", geometric=False)

    my_data = {"density": YTArray(np.random.random(10), "g/cm**3")}
    fake_ds = {"current_time": YTQuantity(10, "Myr")}
    fn = save_as_dataset(fake_ds, "random_data.h5", my_data)
    new_ds = load(fn)
    assert isinstance(new_ds, YTNonspatialDataset)
    yield YTDataFieldTest(fn, "density", geometric=False)
    os.chdir(curdir)
    shutil.rmtree(tmpdir)
