"""
Title: test_ahf.py
Purpose: Contains functions for answer testing the AHF frontend
Notes:
    Copyright (c) 2017, yt Development Team.
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.
"""
from collections import OrderedDict

import pytest

from yt.frontends.ahf.api import AHFHalosDataset
from yt.testing import requires_file
from yt.utilities.answer_testing import \
    framework as fw, \
    utils

# Test data
ahf_halos = 'ahf_halos/snap_N64L16_135.parameter'


#============================================
#                   TestAHF
#============================================
@pytest.mark.skipif(not pytest.config.getvalue('--with-answer-testing'),
    reason="--with-answer-testing not set.")
class TestAHF(fw.AnswerTest):
    #-----
    # test_AHFHalosDataset
    #-----
    @requires_file(ahf_halos)
    def test_AHFHalosDataset(self, ds_ahf_halos):
        assert isinstance(ds_ahf_halos, AHFHalosDataset)

    #-----
    # test_fields_ahf_halos
    #-----
    @utils.requires_ds(ahf_halos)
    def test_fields_ahf_halos(self, ds_ahf_halos):
        fields = ('particle_position_x', 'particle_position_y',
                   'particle_position_z', 'particle_mass')
        hashes = OrderedDict()
        hashes['field_values'] = OrderedDict()
        for field in fields:
            fv_hd = utils.generate_hash(
                self.field_values_test(ds_ahf_halos, field, particle_type=True)
            )
            hashes['field_values'][field] = fv_hd
        # Add function name to hashes
        hashes = {'fields_ahf_halos' : hashes}
        # Save or compare answers
        utils.handle_hashes(self.save_dir, self.answer_file, hashes, self.answer_store)
