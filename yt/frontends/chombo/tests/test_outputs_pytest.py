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

from yt.frontends.chombo.api import \
    ChomboDataset, \
    Orion2Dataset, \
    PlutoDataset
from yt.testing import \
    requires_file, \
    units_override_check
from yt.utilities.answer_testing.answer_tests import small_patch_amr
from yt.utilities.answer_testing.utils import requires_ds

# Test data
gc = "GaussianCloud/data.0077.3d.hdf5"
tb = "TurbBoxLowRes/data.0005.3d.hdf5"
iso = "IsothermalSphere/data.0000.3d.hdf5"
zp = "ZeldovichPancake/plt32.2d.hdf5"
kho = "KelvinHelmholtz/data.0004.hdf5"


@pytest.mark.answer_test
@pytest.mark.usefixtures('answer_file')
class TestChombo:

    @pytest.mark.usefixtures('hashing')
    @pytest.mark.parametrize('ds', [gc], indirect=True)
    def test_gc(self, a, d, w, f, ds):
        self.hashes.update(small_patch_amr(ds, f, w, a, d))

    @pytest.mark.usefixtures('hashing')
    @pytest.mark.parametrize('ds', [tb], indirect=True)
    def test_tb(self, a, d, w, f, ds):
        self.hashes.update(small_patch_amr(ds, f, w, a, d))

    @pytest.mark.usefixtures('hashing')
    @pytest.mark.parametrize('ds', [iso], indirect=True)
    def test_iso(self, a, d, w, f, ds):
        self.hashes.update(small_patch_amr(ds, f, w, a, d))

    @pytest.mark.usefixtures('hashing')
    @pytest.mark.parametrize('ds', [zp], indirect=True)
    def test_zp(self, a, d, w, f, ds):
        self.hashes.update(small_patch_amr(ds, f, w, a, d))

    @pytest.mark.usefixtures('hashing')
    @pytest.mark.parametrize('ds', [kho], indirect=True)
    def test_kho(self, a, d, w, f, ds):
        self.hashes.update(small_patch_amr(ds, f, w, a, d))

    @pytest.mark.parametrize('ds', [zp], indirect=True)
    def test_ChomboDataset(self, ds):
        assert isinstance(ds, ChomboDataset)

    @pytest.mark.parametrize('ds', [gc], indirect=True)
    def test_Orion2Dataset(self, ds):
        assert isinstance(ds, Orion2Dataset)

    @pytest.mark.parametrize('ds', [kho], indirect=True)
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
