"""
Title: test_tipsy.py
Purpose: Tipsy tests using the AGORA dataset
Notes:
    Copyright (c) 2013, yt Development Team.
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.
"""
import pytest

from yt.testing import \
    assert_equal, \
    requires_file
from yt.frontends.tipsy.api import TipsyDataset
import yt.utilities.answer_testing.framework as fw
from yt.utilities.answer_testing import utils


# Test data
pkdgrav = "halo1e11_run1.00400/halo1e11_run1.00400"
gasoline_dmonly = "agora_1e11.00400/agora_1e11.00400"
tipsy_gal = 'TipsyGalaxy/galaxy.00300'




#============================================
#                 TestTipsy
#============================================
@pytest.mark.answer_test
@pytest.mark.usefixtures('answer_file')
class TestTipsy(fw.AnswerTest):
    #-----
    # test_pkdgrav
    #-----
    @pytest.mark.big_data
    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(pkdgrav, file_check = True)
    def test_pkdgrav(self, f, a, d, w, ds_pkdgrav):
        # ds = ds_pkdgrav
        dd = ds_pkdgrav.all_data()
        assert_equal(dd["Coordinates"].shape, (26847360, 3))
        tot = sum(dd[ptype,"Coordinates"].shape[0]
                  for ptype in ds_pkdgrav.particle_types if ptype != "all")
        assert_equal(tot, 26847360)
        ppv_hd = self.pixelized_projection_values_test(ds_pkdgrav, a, f, w, d)
        self.hashes.update({'pixelized_projection_values' : ppv_hd})
        fv_hd = self.field_values_test(ds_pkdgrav, f, d)
        self.hashes.update({'field_values' : fv_hd})
        dobj = utils.create_obj(ds_pkdgrav, d)
        s1 = dobj["ones"].sum()
        s2 = sum(mask.sum() for block, mask in dobj.blocks)
        assert_equal(s1, s2)

    #-----
    # test_gasoline_dmonly
    #-----
    @pytest.mark.big_data
    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(gasoline_dmonly, file_check = True)
    def test_gasoline_dmonly(self, f, a, d, w, ds_gasoline_dmonly):
        ds = ds_gasoline_dmonly
        dd = ds.all_data()
        assert_equal(dd["Coordinates"].shape, (10550576, 3))
        tot = sum(dd[ptype,"Coordinates"].shape[0]
                  for ptype in ds.particle_types if ptype != "all")
        assert_equal(tot, 10550576)
        ppv_hd = self.pixelized_projection_values_test(ds, a, f, w, d)
        self.hashes.update({'pixelized_projection_values' : ppv_hd})
        fv_hd = self.field_values_test(ds, f, d)
        self.hashes.update({'field_values' : fv_hd})
        dobj = utils.create_obj(ds, d)
        s1 = dobj["ones"].sum()
        s2 = sum(mask.sum() for block, mask in dobj.blocks)
        assert_equal(s1, s2)

    #-----
    # test_tipsy_galaxy
    #-----
    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(tipsy_gal)
    def test_tipsy_galaxy(self, f, w, d, a, ds_tipsy_gal):
        self.hashes.update(self.sph_answer(ds_tipsy_gal, 'galaxy.00300', 315372, f, w, d, a))

    #-----
    # test_TipsyDataset
    #-----
    @requires_file(gasoline_dmonly)
    @requires_file(pkdgrav)
    def test_TipsyDataset(self, ds_gasoline_dmonly, ds_pkdgrav):
        assert isinstance(ds_pkdgrav, TipsyDataset)
        assert isinstance(ds_gasoline_dmonly, TipsyDataset)
