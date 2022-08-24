import numpy as np

from yt.frontends.halo_catalog.data_structures import YTHaloCatalogDataset
from yt.frontends.ytdata.utilities import save_as_dataset
from yt.loaders import load as yt_load
from yt.testing import (
    TempDirTest,
    assert_array_equal,
    assert_equal,
    requires_file,
    requires_module,
)
from yt.units.yt_array import YTArray, YTQuantity
from yt.utilities.answer_testing.framework import data_dir_load


def fake_halo_catalog(data):
    filename = "catalog.0.h5"

    ftypes = {field: "." for field in data}
    extra_attrs = {"data_type": "halo_catalog", "num_halos": data["particle_mass"].size}

    ds = {
        "cosmological_simulation": 1,
        "omega_lambda": 0.7,
        "omega_matter": 0.3,
        "hubble_constant": 0.7,
        "current_redshift": 0,
        "current_time": YTQuantity(1, "yr"),
        "domain_left_edge": YTArray(np.zeros(3), "cm"),
        "domain_right_edge": YTArray(np.ones(3), "cm"),
    }
    save_as_dataset(ds, filename, data, field_types=ftypes, extra_attrs=extra_attrs)
    return filename


class HaloCatalogTest(TempDirTest):
    @requires_module("h5py")
    def test_halo_catalog(self):
        rs = np.random.RandomState(3670474)
        n_halos = 100
        fields = ["particle_mass"] + [f"particle_position_{ax}" for ax in "xyz"]
        units = ["g"] + ["cm"] * 3
        data = {
            field: YTArray(rs.random_sample(n_halos), unit)
            for field, unit in zip(fields, units)
        }

        fn = fake_halo_catalog(data)
        ds = yt_load(fn)

        assert type(ds) is YTHaloCatalogDataset

        for field in fields:
            f1 = data[field].in_base()
            f1.sort()
            f2 = ds.r[("all", field)].in_base()
            f2.sort()
            assert_array_equal(f1, f2)

    @requires_module("h5py")
    def test_halo_catalog_boundary_particles(self):
        rs = np.random.RandomState(3670474)
        n_halos = 100
        fields = ["particle_mass"] + [f"particle_position_{ax}" for ax in "xyz"]
        units = ["g"] + ["cm"] * 3
        data = {
            field: YTArray(rs.random_sample(n_halos), unit)
            for field, unit in zip(fields, units)
        }

        data["particle_position_x"][0] = 1.0
        data["particle_position_x"][1] = 0.0
        data["particle_position_y"][2] = 1.0
        data["particle_position_y"][3] = 0.0
        data["particle_position_z"][4] = 1.0
        data["particle_position_z"][5] = 0.0

        fn = fake_halo_catalog(data)
        ds = yt_load(fn)

        assert type(ds) is YTHaloCatalogDataset

        for field in fields:
            f1 = data[field].in_base()
            f1.sort()
            f2 = ds.r[("all", field)].in_base()
            f2.sort()
            assert_array_equal(f1, f2)


t46 = "tiny_fof_halos/DD0046/DD0046.0.h5"


@requires_file(t46)
@requires_module("h5py")
def test_halo_quantities():
    ds = data_dir_load(t46)
    ad = ds.all_data()
    for i in range(ds.index.total_particles):
        hid = int(ad["halos", "particle_identifier"][i])
        halo = ds.halo("halos", hid)
        for field in ["mass", "position", "velocity"]:
            v1 = ad["halos", f"particle_{field}"][i]
            v2 = getattr(halo, field)
            assert_equal(v1, v2, err_msg=f"Halo {hid} {field} field mismatch.")


@requires_file(t46)
@requires_module("h5py")
def test_halo_particles():
    ds = data_dir_load(t46)
    i = ds.r["halos", "particle_mass"].argmax()
    hid = int(ds.r["halos", "particle_identifier"][i])
    halo = ds.halo("halos", hid)
    ids = halo["halos", "member_ids"]
    assert_equal(ids.size, 420)
    assert_equal(ids.min(), 19478.0)
    assert_equal(ids.max(), 31669.0)
