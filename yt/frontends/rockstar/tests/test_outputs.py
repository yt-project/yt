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
    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(r1)
    def test_fields_r1(self, f, ds_r1):
        fv_hd = self.field_values_test(ds_r1, f, particle_type=True)
        self.hashes.update({'field_values' : fv_hd})

    #-----
    # test_RockstarDataset
    #-----
    @requires_file(r1)
    def test_RockstarDataset(self, ds_r1):
        assert isinstance(ds_r1, RockstarDataset)
