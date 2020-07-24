"""
Title: test_fits.py
Purpose: FITS frontend tests
Notes:
    Copyright (c) 2013, yt Development Team.
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.
"""
import pytest

from yt.frontends.fits.data_structures import (
    FITSDataset,
    SpectralCubeFITSDataset,
    SkyDataFITSDataset,
    EventsFITSDataset,
)
from yt.testing import requires_file, units_override_check
from yt.utilities.answer_testing.answer_tests import small_patch_amr
from yt.utilities.answer_testing.utils import requires_ds


# Test data
grs = "radio_fits/grs-50-cube.fits"
vf = "UnigridData/velocity_field_20.fits"
acis = "xray_fits/acisf05356N003_evt2.fits"
A2052 = "xray_fits/A2052_merged_0.3-2_match-core_tmap_bgecorr.fits"


@pytest.mark.answer_test
@pytest.mark.usefixtures("answer_file")
class TestFits:
    @pytest.mark.usefixtures("hashing")
    @requires_ds(grs)
    def test_grs(self, f, a, d, w, ds_grs):
        self.hashes.update(small_patch_amr(ds_grs, f, w, a, d))

    @pytest.mark.usefixtures("hashing")
    @requires_ds(vf)
    def test_velocity_field(self, f, a, d, w, ds_vf):
        self.hashes.update(small_patch_amr(ds_vf, f, w, a, d))

    @pytest.mark.usefixtures("hashing")
    @requires_ds(acis)
    def test_acis(self, f, a, d, w, ds_acis):
        self.hashes.update(small_patch_amr(ds_acis, f, w, a, d))

    @pytest.mark.usefixtures("hashing")
    @requires_ds(A2052)
    def test_A2052(self, f, a, d, w, ds_A2052):
        self.hashes.update(small_patch_amr(ds_A2052, f, w, a, d))

    @requires_file(vf)
    def test_units_override(self):
        units_override_check(vf)

    @requires_file(vf)
    def test_FITSDataset(self, ds_vf):
        assert isinstance(ds_vf, FITSDataset)

    @requires_file(grs)
    def test_SpectralCubeFITSDataset(self, ds_grs):
        assert isinstance(ds_grs, SpectralCubeFITSDataset)

    @requires_file(acis)
    def test_EventsFITSDataset(self, ds_acis):
        assert isinstance(ds_acis, EventsFITSDataset)

    @requires_file(A2052)
    def test_SkyDataFITSDataset(self, ds_A2052):
        assert isinstance(ds_A2052, SkyDataFITSDataset)
