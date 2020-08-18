"""
title: test_chombo.py
Purpose: Chombo frontend tests
Notes:
    Copyright (c) 2013, yt Development Team.
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.
"""
import pytest

from yt.frontends.chombo.api import ChomboDataset, Orion2Dataset, PlutoDataset
from yt.testing import requires_file, units_override_check
from yt.utilities.answer_testing.answer_tests import small_patch_amr

# Test data
gc = "GaussianCloud/data.0077.3d.hdf5"
tb = "TurbBoxLowRes/data.0005.3d.hdf5"
iso = "IsothermalSphere/data.0000.3d.hdf5"
zp = "ZeldovichPancake/plt32.2d.hdf5"
kho = "KelvinHelmholtz/data.0004.hdf5"


ds_list = [
    gc,
    tb,
    iso,
    zp,
    kho,
]
a_list = [0, 1, 2]
w_list = [None, "density"]
w_zp = [None, "rhs"]
d_list = [None, ("sphere", ("max", (0.1, "unitary")))]
d_zp = [(None, ("sphere", ("c", (0.1, "unitary")))]
f_list = ["density", "velocity_magnitude", "magnetic_field_x"]
f_zp = ["rhs", "phi"]

pairs_list = [
    [gc, f_list, d_list, w_list],
    [tb, f_list, d_list, w_list],
    [iso, f_list, d_list, w_list],
    [zp, f_zp, d_zp, w_zp],
    [kho, f_list, d_list, w_list],

gv_pairs = [
    (i[0], f) for i in ds_list for f in i[1]] 
]
fv_pairs = [
    (i[0], f, d) for i in pair_list for f in i[1] for d in i[2]]
]
pv_pairs = [
    (i[0], f, d, w) for i in pair_list for f in i[1] for d in i[2] for w in i[3]]
]


@pytest.mark.answer_test
class TestChombo:
    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds", ds_list, indirect=True)
    def test_gh_pr(self, ds):
        self.hashes.update({"grid_hierarchy": grid_hierarchy(ds)})
        self.hashes.update({"parentage_relationships": parentage_relationships(ds)})

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds, f", gv_pairs, indirect=True)
    def test_gv(self, f, ds):
        self.hashes.update({"grid_values": grid_values(ds, f)})

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds, f, d", fv_pairs, indirect=True)
    def test_fv(self, d, f, ds):
        self.hashes.update({"field_values": field_values(ds, f, d)})

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds, f, d, w", pv_pairs, indirect=True)
    @pytest.mark.parametrize("a", a_list, indirect=True)
    def test_pv(self, a, d, w, f, ds):
        self.hashes.update({"projection_values": projection_values(ds, a, f, w, d)})

    @pytest.mark.parametrize("ds", [zp], indirect=True)
    def test_ChomboDataset(self, ds):
        assert isinstance(ds, ChomboDataset)

    @pytest.mark.parametrize("ds", [gc], indirect=True)
    def test_Orion2Dataset(self, ds):
        assert isinstance(ds, Orion2Dataset)

    @pytest.mark.parametrize("ds", [kho], indirect=True)
    def test_PlutoDataset(self, ds):
        assert isinstance(ds, PlutoDataset)

    @requires_file(zp)
    def test_units_override_zp(self):
        units_override_check(zp)

    @requires_file(gc)
    def test_units_override_gc(self):
        units_override_check(gc)

    @requires_file(kho)
    def test_units_override_kho(self):
        units_override_check(kho)
