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
    EventsFITSDataset,
    FITSDataset,
    SkyDataFITSDataset,
    SpectralCubeFITSDataset,
)
from yt.frontends.fits.misc import setup_counts_fields
from yt.testing import requires_file, units_override_check
from yt.utilities.answer_testing.answer_tests import (
    field_values,
    grid_hierarchy,
    grid_values,
    parentage_relationships,
    projection_values,
)
from yt.utilities.answer_testing.utils import requires_ds
from yt.utilities.exceptions import YTOutputNotIdentified


# Test data
grs = "radio_fits/grs-50-cube.fits"
vf = "UnigridData/velocity_field_20.fits"
acis = "xray_fits/acisf05356N003_evt2.fits"
A2052 = "xray_fits/A2052_merged_0.3-2_match-core_tmap_bgecorr.fits"


def get_acis():
    try:
        ds = data_dir_load(acis, cls=EventsFITSDataset)
        ebounds = [(0.1, 2.0), (2.0, 5.0)]
        setup_counts_fields(ds, ebounds)
        return ds
    except YTOutputNotIdentified:
        return pytest.skip("acis not found.")


grs_kwargs = {"kwargs" : {"nan_mask" : 0.0}, "cls" : SpectralCubeFITSDataset}
vf_kwargs = {"cls" : FITSDataset}
A2052_kwargs = {"cls" : SkyDataFITSDataset}

_fields_grs = ("temperature",)
_fields_vels = ("velocity_x", "velocity_y", "velocity_z")
_fields_acis = ("counts_0.1-2.0", "counts_2.0-5.0")
_fields_A2052 = ("flux",)

a_list = [0, 1, 2]
d_list = [None, ("sphere", ("c", (0.1, "unitary")))]
w_list = [None, "ones"]
f_list = [
    _fields_grs,
    _fields_vels,
    _fields_acis,
    _fields_A2052,
]
ds_list = [
    [grs, grs_kwargs], 
    [vf, vf_kwargs],
    get_acis(),
    [A2052, A2052_kwargs],
]


def get_pairs():
    pairs = []
    for i, ds in enumerate(ds_list):
        for f in f_list[i]:
            pairs.append((ds, f))
    return pairs


@pytest.mark.answer_test
class TestFits:
    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds", ds_list, indirect=True)
    def test_gh_pr(self, ds):
        self.hashes.update({"grid_hierarchy": grid_hierarchy(ds)})
        self.hashes.update({"parentage_relationships": parentage_relationships(ds)})

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds, f", get_pairs(), indirect=True)
    def test_gv(self, f, ds):
        self.hashes.update({"grid_values": grid_values(ds, f)})

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds, f", get_pairs(), indirect=True)
    @pytest.mark.parametrize("d", d_list, indirect=True)
    def test_fv(self, d, f, ds):
        self.hashes.update({"field_values": field_values(ds, f, d)})

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds, f", get_pairs(), indirect=True)
    @pytest.mark.parametrize("d", d_list, indirect=True)
    @pytest.mark.parametrize("a", a_list, indirect=True)
    @pytest.mark.parametrize("w", w_list, indirect=True)
    def test_pv(self, a, d, w, f, ds):
        self.hashes.update({"projection_values": projection_values(ds, a, f, w, d)})

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
