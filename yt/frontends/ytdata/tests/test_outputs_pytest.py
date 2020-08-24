import os

import numpy as np
import pytest

from yt.convenience import load
from yt.data_objects.api import create_profile
from yt.frontends.ytdata.api import (
    YTDataContainerDataset,
    YTGridDataset,
    YTNonspatialDataset,
    YTProfileDataset,
    YTSpatialPlotDataset,
    save_as_dataset,
)
from yt.testing import assert_array_equal
from yt.units.yt_array import YTArray, YTQuantity
from yt.utilities.answer_testing import utils
from yt.utilities.answer_testing.answer_tests import yt_data_field
from yt.visualization.profile_plotter import PhasePlot, ProfilePlot

# Test data
enzotiny = "enzo_tiny_cosmology/DD0046/DD0046"


@pytest.mark.answer_test
@pytest.mark.usefixtures("temp_dir")
class TestYTData:
    answer_file = None

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds", [enzotiny], indirect=True)
    def test_datacontainer_data(self, field, ds):
        sphere = ds.sphere(ds.domain_center, (10, "Mpc"))
        fn = sphere.save_as_dataset(fields=["density", "particle_mass"])
        full_fn = os.path.join(os.getcwd(), fn)
        sphere_ds = load(full_fn)
        utils.compare_unit_attributes(ds, sphere_ds)
        assert isinstance(sphere_ds, YTDataContainerDataset)
        ytft = yt_data_field(sphere_ds, field, True)
        self.hashes.update({"yt_data_field": ytft})
        cr = ds.cut_region(sphere, ['obj["temperature"] > 1e4'])
        fn = cr.save_as_dataset(fields=["temperature"])
        full_fn = os.path.join(os.getcwd(), fn)
        cr_ds = load(full_fn)
        assert isinstance(cr_ds, YTDataContainerDataset)
        assert (cr["temperature"] == cr_ds.data["temperature"]).all()

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds", [enzotiny], indirect=True)
    def test_covering_grid_datacontainer_data(self, field, ds):
        cg = ds.covering_grid(level=0, left_edge=[0.25] * 3, dims=[16] * 3)
        fn = cg.save_as_dataset(
            fields=["density", "particle_mass", "particle_position"]
        )
        full_fn = os.path.join(os.getcwd(), fn)
        cg_ds = load(full_fn)
        utils.compare_unit_attributes(ds, cg_ds)
        assert isinstance(cg_ds, YTGridDataset)
        assert (
            cg["all", "particle_position"].shape
            == cg_ds.r["all", "particle_position"].shape
        )
        ytft = yt_data_field(cg_ds, field, True)
        self.hashes.update({"yt_data_field": ytft})

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds", [enzotiny], indirect=True)
    def test_arbitrary_grid_datacontainer_data(self, field, ds):
        ag = ds.arbitrary_grid(
            left_edge=[0.25] * 3, right_edge=[0.75] * 3, dims=[16] * 3
        )
        fn = ag.save_as_dataset(fields=["density", "particle_mass"])
        full_fn = os.path.join(os.getcwd(), fn)
        ag_ds = load(full_fn)
        utils.compare_unit_attributes(ds, ag_ds)
        assert isinstance(ag_ds, YTGridDataset)
        ytft = yt_data_field(ag_ds, field, True)
        self.hashes.update({"yt_data_field": ytft})

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds", [enzotiny], indirect=True)
    def test_frb_datacontainer_data(self, field, ds):
        my_proj = ds.proj("density", "x", weight_field="density")
        frb = my_proj.to_frb(1.0, (800, 800))
        fn = frb.save_as_dataset(fields=["density"])
        frb_ds = load(fn)
        assert_array_equal(frb["density"], frb_ds.data["density"])
        utils.compare_unit_attributes(ds, frb_ds)
        assert isinstance(frb_ds, YTGridDataset)
        ytft = yt_data_field(frb_ds, field, False)
        self.hashes.update({"yt_data_field": ytft})

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds", [enzotiny], indirect=True)
    def test_spatial_data(self, field, ds):
        proj = ds.proj("density", "x", weight_field="density")
        fn = proj.save_as_dataset()
        full_fn = os.path.join(os.getcwd(), fn)
        proj_ds = load(full_fn)
        utils.compare_unit_attributes(ds, proj_ds)
        assert isinstance(proj_ds, YTSpatialPlotDataset)
        ytft = yt_data_field(proj_ds, field, False)
        self.hashes.update({"yt_data_field": ytft})

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds", [enzotiny], indirect=True)
    def test_profile_data1(self, field, ds):
        ad = ds.all_data()
        profile_1d = create_profile(
            ad, "density", "temperature", weight_field="cell_mass"
        )
        fn = profile_1d.save_as_dataset()
        full_fn = os.path.join(os.getcwd(), fn)
        prof_1d_ds = load(full_fn)
        utils.compare_unit_attributes(ds, prof_1d_ds)
        assert isinstance(prof_1d_ds, YTProfileDataset)
        for field in profile_1d.standard_deviation:
            assert_array_equal(
                profile_1d.standard_deviation[field],
                prof_1d_ds.profile.standard_deviation["data", field[1]],
            )
        p1 = ProfilePlot(
            prof_1d_ds.data, "density", "temperature", weight_field="cell_mass"
        )
        p1.save()
        ytft = yt_data_field(prof_1d_ds, field, False)
        self.hashes.update({"yt_data_field": ytft})

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds", [enzotiny], indirect=True)
    def test_profile_data2(self, field, ds):
        ad = ds.all_data()
        profile_2d = create_profile(
            ad,
            ["density", "temperature"],
            "cell_mass",
            weight_field=None,
            n_bins=(128, 128),
        )
        fn = profile_2d.save_as_dataset()
        full_fn = os.path.join(os.getcwd(), fn)
        prof_2d_ds = load(full_fn)
        utils.compare_unit_attributes(ds, prof_2d_ds)
        assert isinstance(prof_2d_ds, YTProfileDataset)
        p2 = PhasePlot(
            prof_2d_ds.data, "density", "temperature", "cell_mass", weight_field=None
        )
        p2.save()
        ytft = yt_data_field(prof_2d_ds, field, False)
        self.hashes.update({"yt_data_field": ytft})

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds", [enzotiny], indirect=True)
    def test_nonspatial_data1(self, field, ds):
        region = ds.box([0.25] * 3, [0.75] * 3)
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
        ytft = yt_data_field(array_ds, field, False)
        self.hashes.update({"yt_data_field": ytft})

    @pytest.mark.usefixtures("hashing")
    def test_nonspatial_data2(self, field):
        my_data = {"density": YTArray(np.linspace(1.0, 20.0, 10), "g/cm**3")}
        fake_ds = {"current_time": YTQuantity(10, "Myr")}
        fn = "random_data.h5"
        save_as_dataset(fake_ds, fn, my_data)
        full_fn = os.path.join(os.getcwd(), fn)
        new_ds = load(full_fn)
        assert isinstance(new_ds, YTNonspatialDataset)
        ytft = yt_data_field(new_ds, field, False)
        self.hashes.update({"yt_data_field": ytft})
