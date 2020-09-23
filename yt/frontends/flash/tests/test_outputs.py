import numpy as np
import pytest

from yt.frontends.flash.api import FLASHDataset, FLASHParticleDataset
from yt.testing import (
    ParticleSelectionComparison,
    assert_equal,
    requires_file,
    units_override_check,
)
from yt.utilities.answer_testing.answer_tests import (
    field_values,
    grid_hierarchy,
    grid_values,
    nbody_answer,
    parentage_relationships,
    projection_values,
)

# Test data
sloshing = "GasSloshingLowRes/sloshing_low_res_hdf5_plt_cnt_0300"
wt = "WindTunnel/windtunnel_4lev_hdf5_plt_cnt_0030"
fid_1to3_b1 = "fiducial_1to3_b1/fiducial_1to3_b1_hdf5_part_0080"
dens_turb_mag = "DensTurbMag/DensTurbMag_hdf5_plt_cnt_0015"

a_list = [0, 1, 2]
w_list = [None, "density"]
d_list = [None, ("sphere", ("max", (0.1, "unitary")))]
f_list = [
    ["temperature", "density", "velocity_magnitude"],
    ["temperature", "density"],
]
ds_list = [
    sloshing,
    wt,
]


def get_pairs():
    pairs = []
    for i, d in enumerate(ds_list):
        for f in f_list[i]:
            pairs.append((d, f))
    return pairs


@pytest.mark.answer_test
class TestFlash:
    answer_file = None
    saved_hashes = None

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds", ds_list, indirect=True)
    def test_gh_pr(self, ds, big_data):
        if str(ds) == "sloshing_low_res_hdf5_plt_cnt_0300" and not big_data:
            pytest.skip("--answer-big-data not set.")
        self.hashes.update({"grid_hierarchy": grid_hierarchy(ds)})
        self.hashes.update({"parentage_relationships": parentage_relationships(ds)})

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds, f", get_pairs(), indirect=True)
    def test_gv(self, f, ds, big_data):
        if str(ds) == "sloshing_low_res_hdf5_plt_cnt_0300" and not big_data:
            pytest.skip("--answer-big-data not set.")
        self.hashes.update({"grid_values": grid_values(ds, f)})

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds, f", get_pairs(), indirect=True)
    @pytest.mark.parametrize("d", d_list, indirect=True)
    def test_fv(self, d, f, ds, big_data):
        if str(ds) == "sloshing_low_res_hdf5_plt_cnt_0300" and not big_data:
            pytest.skip("--answer-big-data not set.")
        self.hashes.update({"field_values": field_values(ds, f, d)})

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds, f", get_pairs(), indirect=True)
    @pytest.mark.parametrize("d", d_list, indirect=True)
    @pytest.mark.parametrize("a", a_list, indirect=True)
    @pytest.mark.parametrize("w", w_list, indirect=True)
    def test_pv(self, a, d, w, f, ds, big_data):
        if str(ds) == "sloshing_low_res_hdf5_plt_cnt_0300" and not big_data:
            pytest.skip("--answer-big-data not set.")
        self.hashes.update({"projection_values": projection_values(ds, a, f, w, d)})

    @pytest.mark.big_data
    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds", [fid_1to3_b1], indirect=True)
    def test_fid_1to3_b1(self, f, w, d, a, ds):
        self.hashes.update(
            nbody_answer(ds, "fiducial_1to3_b1_hdf5_part_0080", 6684119, f, w, d, a)
        )

    @pytest.mark.parametrize("ds", [wt], indirect=True)
    def test_FLASHDataset(self, ds):
        assert isinstance(ds, FLASHDataset)

    @requires_file(sloshing)
    def test_units_override(self):
        units_override_check(sloshing)

    @pytest.mark.parametrize("ds", [fid_1to3_b1], indirect=True)
    def test_FLASHParticleDataset(self, ds):
        assert isinstance(ds, FLASHParticleDataset)

    @pytest.mark.parametrize("ds", [dens_turb_mag], indirect=True)
    def test_FLASH25_dataset(self, ds):
        assert_equal(ds.parameters["time"], 751000000000.0)
        assert_equal(ds.domain_dimensions, np.array([8, 8, 8]))
        assert_equal(ds.domain_left_edge, ds.arr([-2e18, -2e18, -2e18], "code_length"))
        assert_equal(ds.index.num_grids, 73)
        dd = ds.all_data()
        dd["density"]

    @pytest.mark.parametrize("ds", [fid_1to3_b1], indirect=True)
    def test_FLASHParticleDataset_selection(self, ds):
        psc = ParticleSelectionComparison(ds)
        psc.run_defaults()

    @pytest.mark.parametrize("ds", [sloshing], indirect=True)
    def test_mu(self, ds):
        sp = ds.sphere("c", (0.1, "unitary"))
        assert np.all(
            sp["gas", "mean_molecular_weight"] == ds.parameters["eos_singlespeciesa"]
        )
