"""
Title: test_artio.py
Purpose: Contains ARTIO frontend tests
Notes:
    Copyright (c) 2013, yt Development Team.
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.
"""
import pytest

from yt.frontends.artio.api import ARTIODataset
from yt.testing import \
    assert_allclose_units, \
    assert_equal, \
    requires_file, \
    units_override_check
import yt.utilities.answer_testing.framework as fw
from yt.utilities.answer_testing import utils


# Data file
sizmbhloz = "sizmbhloz-clref04SNth-rs9_a0.9011/"
sizmbhloz += "sizmbhloz-clref04SNth-rs9_a0.9011.art"


#============================================
#                 TestArtIo
#============================================
@pytest.mark.answer_test
@pytest.mark.usefixtures('answer_file')
class TestArtIo(fw.AnswerTest):
    #-----
    # test_sizmbhloz
    #-----
    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(sizmbhloz)
    def test_sizmbhloz(self, d, a, w, f, ds_sizmbhloz):
        ds_sizmbhloz.max_range = 1024*1024
        # Run tests
        fv_hd = self.field_values_test(ds_sizmbhloz, f, d)
        self.hashes.update({'field_values_test' : fv_hd}) 
        ppv_hd = self.pixelized_projection_values_test(ds_sizmbhloz, a, f, w, d)
        self.hashes.update({'pixelized_projection_values' : ppv_hd}) 
        dobj = utils.create_obj(ds_sizmbhloz, d)
        s1 = dobj["ones"].sum()
        s2 = sum(mask.sum() for block, mask in dobj.blocks)
        assert_equal(s1, s2)
        assert_equal(ds_sizmbhloz.particle_type_counts, {'N-BODY': 100000, 'STAR': 110650})

    #-----
    # test_ARTIODataset
    #-----
    @requires_file(sizmbhloz)
    def test_ARTIODataset(self, ds_sizmbhloz):
        assert isinstance(ds_sizmbhloz, ARTIODataset)

    #-----
    # test_units_override
    #-----
    @requires_file(sizmbhloz)
    def test_units_override(self, ds_sizmbhloz):
        units_override_check(ds_sizmbhloz, sizmbhloz)

    #-----
    # test_particle_derived_field
    #-----
    @requires_file(sizmbhloz)
    def test_particle_derived_field(self, ds_sizmbhloz):
        def star_age_alias(field, data):
            # test to make sure we get back data in the correct units
            # during field detection
            return data['STAR', 'age'].in_units('Myr')

        ds_sizmbhloz.add_field(("STAR", "new_field"), function=star_age_alias,
                     units='Myr', sampling_type="particle")
        ad = ds_sizmbhloz.all_data()
        assert_allclose_units(ad['STAR', 'age'].in_units("Myr"),
                              ad["STAR", "new_field"])
