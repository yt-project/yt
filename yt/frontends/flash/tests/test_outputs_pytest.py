"""
Title: test_flash.py
Purpose: FLASH frontend tests
Notes:
    Copyright (c) 2013, yt Development Team.
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.
"""
import numpy as np
import pytest

from yt.frontends.flash.api import FLASHDataset, FLASHParticleDataset
from yt.testing import (
    assert_equal,
    ParticleSelectionComparison,
    requires_file,
    units_override_check,
)
from yt.utilities.answer_testing.answer_tests import nbody_answer, small_patch_amr

# Test data
sloshing = "GasSloshingLowRes/sloshing_low_res_hdf5_plt_cnt_0300"
wt = "WindTunnel/windtunnel_4lev_hdf5_plt_cnt_0030"
fid_1to3_b1 = "fiducial_1to3_b1/fiducial_1to3_b1_hdf5_part_0080"
dens_turb_mag = "DensTurbMag/DensTurbMag_hdf5_plt_cnt_0015"


@pytest.mark.answer_test
@pytest.mark.usefixtures("answer_file")
class TestFlash:
    @pytest.mark.big_data
    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds", [sloshing], indirect=True)
    def test_sloshing(self, a, d, w, f, ds):
        self.hashes.update(small_patch_amr(ds, f, w, a, d))

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds", [wt], indirect=True)
    def test_wind_tunnel(self, a, d, w, f, ds):
        self.hashes.update(small_patch_amr(ds, f, w, a, d))

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

    @pytest.mark.big_data
    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds", [fid_1to3_b1], indirect=True)
    def test_fid_1to3_b1(self, f, w, d, a, ds):
        self.hashes.update(
            nbody_answer(ds, "fiducial_1to3_b1_hdf5_part_0080", 6684119, f, w, d, a)
        )

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
