"""
Title: test_rockstar.py
Purpose: Rockstar frontend tests using rockstar_halos dataset
Notes:
    Copyright (c) 2015, yt Development Team.
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.
"""
from collections import OrderedDict

import pytest

from yt.testing import requires_file
from yt.frontends.rockstar.api import RockstarDataset
import yt.utilities.answer_testing.framework as fw
from yt.utilities.answer_testing import utils


# Test data
r1 = "rockstar_halos/halos_0.0.bin"


# Globals
_fields = ("particle_position_x", "particle_position_y",
           "particle_position_z", "particle_mass")


#============================================
#               TestRockstar
#============================================
@pytest.mark.skipif(not pytest.config.getvalue('--with-answer-testing'),
    reason="--with-answer-testing not set.")
@pytest.mark.usefixtures('answer_file')
class TestRockstar(fw.AnswerTest):
    #-----
    # test_fields_r1
    #-----
    @utils.requires_ds(r1)
    def test_fields_r1(self, ds_r1):
        hashes = OrderedDict()
        hashes['field_values'] = OrderedDict()
        for field in _fields:
            fv_hd = utils.generate_hash(
                self.field_values_test(ds_r1, field, particle_type=True)
            )
            hashes['field_values'][field] = fv_hd
        hashes = {'fields_r1' : hashes}
        utils.handle_hashes(self.save_dir, self.answer_file, hashes, self.answer_store)

    #-----
    # test_RockstarDataset
    #-----
    @requires_file(r1)
    def test_RockstarDataset(self, ds_r1):
        assert isinstance(ds_r1, RockstarDataset)
