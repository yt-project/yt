"""
title: test_ytdata.py
Purpose: ytdata frontend tests using enzo_tiny_cosmology
Notes:
    Copyright (c) yt Development Team. All rights reserved.
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.
"""
import os
from collections import OrderedDict

import numpy as np
import pytest

from yt.convenience import \
    load
from yt.data_objects.api import \
    create_profile
from yt.frontends.ytdata.api import \
    YTDataContainerDataset, \
    YTSpatialPlotDataset, \
    YTGridDataset, \
    YTNonspatialDataset, \
    YTProfileDataset, \
    save_as_dataset
from yt.testing import \
    assert_array_equal, \
    assert_fname, \
    fake_random_ds, \
    requires_module
from yt.units.yt_array import \
    YTArray, \
    YTQuantity
from yt.visualization.plot_window import \
    SlicePlot, \
    ProjectionPlot
from yt.visualization.profile_plotter import \
    ProfilePlot, \
    PhasePlot
import yt.utilities.answer_testing.framework as fw
from yt.utilities.answer_testing import utils

# Test data
enzotiny = "enzo_tiny_cosmology/DD0046/DD0046"


#============================================
#                TestYTData
#============================================
@pytest.mark.skipif(not pytest.config.getvalue('--with-answer-testing'),
    reason="--with-answer-testing not set.")
@pytest.mark.usefixtures('temp_dir', 'answer_file')
class TestYTData(fw.AnswerTest):
    #-----
    # test_datacontainer_data
    #-----
    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(enzotiny)
    def test_datacontainer_data(self, field, ds_enzotiny):
        ds = ds_enzotiny
        sphere = ds.sphere(ds.domain_center, (10, "Mpc"))
        fn = sphere.save_as_dataset(fields=["density", "particle_mass"])
        full_fn = os.path.join(os.getcwd(), fn)
        sphere_ds = load(full_fn)
        utils.compare_unit_attributes(ds, sphere_ds)
        assert isinstance(sphere_ds, YTDataContainerDataset)
        ytft_hd = self.yt_field_test(sphere_ds, field, True)
        self.hashes.update({'yt_field_test' : ytft_hd})
        cr = ds.cut_region(sphere, ['obj["temperature"] > 1e4'])
        fn = cr.save_as_dataset(fields=["temperature"])
        full_fn = os.path.join(os.getcwd(), fn)
        cr_ds = load(full_fn)
        assert isinstance(cr_ds, YTDataContainerDataset)
        assert (cr["temperature"] == cr_ds.data["temperature"]).all()

    #-----
    # test_covering_grid_datacontainer_data
    #-----
    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(enzotiny)
    def test_covering_grid_datacontainer_data(self, field, ds_enzotiny):
        ds = ds_enzotiny
        cg = ds.covering_grid(level=0, left_edge=[0.25]*3, dims=[16]*3)
        fn = cg.save_as_dataset(fields=["density", "particle_mass"])
        full_fn = os.path.join(os.getcwd(), fn)
        cg_ds = load(full_fn)
        utils.compare_unit_attributes(ds, cg_ds)
        assert isinstance(cg_ds, YTGridDataset)
        ytft_hd = self.yt_field_test(cg_ds, field, True)
        self.hashes.update({'yt_field_test' :  ytft_hd})

    #-----
    # test_arbitrary_grid_datacontainer_data
    #-----
    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(enzotiny)
    def test_arbitrary_grid_datacontainer_data(self, field, ds_enzotiny):
        ds = ds_enzotiny
        ag = ds.arbitrary_grid(left_edge=[0.25]*3, right_edge=[0.75]*3,
                               dims=[16]*3)
        fn = ag.save_as_dataset(fields=["density", "particle_mass"])
        full_fn = os.path.join(os.getcwd(), fn)
        ag_ds = load(full_fn)
        utils.compare_unit_attributes(ds, ag_ds)
        assert isinstance(ag_ds, YTGridDataset)
        ytft_hd = self.yt_field_test(ag_ds, field, True)
        self.hashes.update({'yt_field_test' : ytft_hd})

    #-----
    # test_frb_datacontainer_data
    #-----
    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(enzotiny)
    def test_frb_datacontainer_data(self, field, ds_enzotiny):
        ds = ds_enzotiny
        my_proj = ds.proj("density", "x", weight_field="density")
        frb = my_proj.to_frb(1.0, (800, 800))
        fn = frb.save_as_dataset(fields=["density"])
        frb_ds = load(fn)
        assert_array_equal(frb["density"], frb_ds.data["density"])
        utils.compare_unit_attributes(ds, frb_ds)
        assert isinstance(frb_ds, YTGridDataset)
        ytft_hd = self.yt_field_test(frb_ds, field, False)
        self.hashes.update({'yt_field_test' : ytft_hd})

    #-----
    # test_spatial_data
    #-----
    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(enzotiny)
    def test_spatial_data(self, field, ds_enzotiny):
        ds = ds_enzotiny
        proj = ds.proj("density", "x", weight_field="density")
        fn = proj.save_as_dataset()
        full_fn = os.path.join(os.getcwd(), fn)
        proj_ds = load(full_fn)
        utils.compare_unit_attributes(ds, proj_ds)
        assert isinstance(proj_ds, YTSpatialPlotDataset)
        ytft_hd = self.yt_field_test(proj_ds, field, False)
        self.hashes.update({'yt_field_test' : ytft_hd})

    #-----
    # test_profile_data
    #-----
    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(enzotiny)
    def test_profile_data1(self, field, ds_enzotiny):
        ds = ds_enzotiny
        ad = ds.all_data()
        profile_1d = create_profile(ad, "density", "temperature",
                                    weight_field="cell_mass")
        fn = profile_1d.save_as_dataset()
        full_fn = os.path.join(os.getcwd(), fn)
        prof_1d_ds = load(full_fn)
        utils.compare_unit_attributes(ds, prof_1d_ds)
        assert isinstance(prof_1d_ds, YTProfileDataset)
        for field in profile_1d.standard_deviation:
            assert_array_equal(
                profile_1d.standard_deviation[field],
                prof_1d_ds.profile.standard_deviation['data', field[1]])
        p1 = ProfilePlot(prof_1d_ds.data, "density", "temperature",
                         weight_field="cell_mass")
        p1.save()
        ytft_hd = self.yt_field_test(prof_1d_ds, field, False)
        self.hashes.update({'yt_field_test' : ytft_hd})

    #-----
    # test_profile_data2
    #-----
    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(enzotiny)
    def test_profile_data2(self, field, ds_enzotiny):
        ds = ds_enzotiny
        profile_2d = create_profile(ad, ["density", "temperature"],
                                   "cell_mass", weight_field=None,
                                   n_bins=(128, 128))
        fn = profile_2d.save_as_dataset()
        full_fn = os.path.join(os.getcwd(), fn)
        prof_2d_ds = load(full_fn)
        utils.compare_unit_attributes(ds, prof_2d_ds)
        assert isinstance(prof_2d_ds, YTProfileDataset)
        p2 = PhasePlot(prof_2d_ds.data, "density", "temperature",
                       "cell_mass", weight_field=None)
        p2.save()
        ytft_hd = self.yt_field_test(prof_2d_ds, field, False)
        self.hashes.update({'yt_field_test' : ytft_hd})

    #-----
    # test_nonspatial_data1
    #-----
    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(enzotiny)
    def test_nonspatial_data1(self, field, ds_enzotiny):
        ds = ds_enzotiny
        region = ds.box([0.25]*3, [0.75]*3)
        sphere = ds.sphere(ds.domain_center, (10, "Mpc"))
        my_data = {}
        my_data["region_density"] = region["density"]
        my_data["sphere_density"] = sphere["density"]
        fn = "test_data.h5"
        save_as_dataset(ds, fn, my_data)
        full_fn = os.path.join(os.getcwd(), fn)
        array_ds = load(full_fn)
        utils.compare_unit_attributes(ds, array_ds)
        assert isinstance(array_ds, YTNonspatialDataset)
        ytft_hd = self.yt_field_test(array_ds, field, False)
        self.hashes.update({'yt_field_test' : ytft_hd})
    
    #-----
    # test_nonspatial_data2
    #-----
    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(enzotiny)
    def test_nonspatial_data2(self, field, ds_enzotiny):
        ds = ds_enzotiny
        my_data = {"density": YTArray(np.linspace(1.,20.,10), "g/cm**3")}
        fake_ds = {"current_time": YTQuantity(10, "Myr")}
        fn = "random_data.h5"
        save_as_dataset(fake_ds, fn, my_data)
        full_fn = os.path.join(os.getcwd(), fn)
        new_ds = load(full_fn)
        assert isinstance(new_ds, YTNonspatialDataset)
        ytft_hd = self.yt_field_test(new_ds, field, False)
        self.hashes.update({'yt_field_test' : ytft_hd})

    #-----
    # test_plot_data
    #-----
    @requires_module('h5py')
    def test_plot_data(self):
        ds = fake_random_ds(16)
        plot = SlicePlot(ds, 'z', 'density')
        fn = plot.data_source.save_as_dataset('slice.h5')
        ds_slice = load(fn)
        p = SlicePlot(ds_slice, 'z', 'density')
        fn = p.save()
        assert_fname(fn[0])
        plot = ProjectionPlot(ds, 'z', 'density')
        fn = plot.data_source.save_as_dataset('proj.h5')
        ds_proj = load(fn)
        p = ProjectionPlot(ds_proj, 'z', 'density')
        fn = p.save()
        assert_fname(fn[0])
        plot = SlicePlot(ds, [1, 1, 1], 'density')
        fn = plot.data_source.save_as_dataset('oas.h5')
        ds_oas = load(fn)
        p = SlicePlot(ds_oas, [1, 1, 1], 'density')
        fn = p.save()
        assert_fname(fn[0])
