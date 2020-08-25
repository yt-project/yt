import numpy as np
import pytest

from yt.convenience import load
from yt.frontends.halo_catalog.data_structures import HaloCatalogDataset
from yt.testing import assert_array_equal, requires_module
from yt.units.yt_array import YTArray
from yt.utilities.answer_testing.utils import fake_halo_catalog


@pytest.mark.answer_test
@pytest.mark.usefixtures("temp_dir")
class TestHaloCatalog:
    answer_file = None

    @requires_module("h5py")
    def test_halo_catalog(self):
        rs = np.random.RandomState(3670474)
        n_halos = 100
        fields = [
            f"particle_{name}" for name in ["mass"] + [f"position_{ax}" for ax in "xyz"]
        ]
        units = ["g"] + ["cm"] * 3
        data = dict(
            (field, YTArray(rs.random_sample(n_halos), unit))
            for field, unit in zip(fields, units)
        )
        fn = fake_halo_catalog(data)
        ds = load(fn)
        assert isinstance(ds, HaloCatalogDataset)
        for field in fields:
            f1 = data[field].in_base()
            f1.sort()
            f2 = ds.r[field].in_base()
            f2.sort()
            assert_array_equal(f1, f2)

    @requires_module("h5py")
    def test_halo_catalog_boundary_particles(self):
        rs = np.random.RandomState(3670474)
        n_halos = 100
        fields = [
            f"particle_{name}" for name in ["mass"] + [f"position_{ax}" for ax in "xyz"]
        ]
        units = ["g"] + ["cm"] * 3
        data = dict(
            (field, YTArray(rs.random_sample(n_halos), unit))
            for field, unit in zip(fields, units)
        )
        data["particle_position_x"][0] = 1.0
        data["particle_position_x"][1] = 0.0
        data["particle_position_y"][2] = 1.0
        data["particle_position_y"][3] = 0.0
        data["particle_position_z"][4] = 1.0
        data["particle_position_z"][5] = 0.0
        fn = fake_halo_catalog(data)
        ds = load(fn)
        assert isinstance(ds, HaloCatalogDataset)
        for field in ["particle_mass"]:
            f1 = data[field].in_base()
            f1.sort()
            f2 = ds.r[field].in_base()
            f2.sort()
            assert_array_equal(f1, f2)
