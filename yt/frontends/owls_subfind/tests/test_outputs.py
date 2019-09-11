"""
Title: test_owls_subfind.py
Purpose: OWLSSubfind frontend tests using owls_fof_halos datasets
Notes:
    Copyright (c) 2013, yt Development Team.
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.
"""
import os.path
from collections import OrderedDict

import pytest

from yt.frontends.owls_subfind.api import OWLSSubfindDataset
from yt.testing import \
    assert_equal, \
    requires_file
import yt.utilities.answer_testing.framework as fw
from yt.utilities.answer_testing import utils

# Test data
g1 = "owls_fof_halos/groups_001/group_001.0.hdf5"
g8 = "owls_fof_halos/groups_008/group_008.0.hdf5"


# Answer file
answer_file = 'owls_subfind_answers.yaml'


#============================================
#              TestOwlsSubfind
#============================================
@pytest.mark.skipif(not pytest.config.getvalue('--with-answer-testing'),
    reason="--with-answer-testing not set.")
class TestOwlsSubfind(fw.AnswerTest):
    #-----
    # test_fields_g8
    #-----
    @utils.requires_ds(g8)
    def test_fields_g8(self, ds_g8):
        hashes = OrderedDict()
        hashes['field_values'] = OrderedDict()
        fields = ("particle_position_x", "particle_position_y",
                   "particle_position_z", "particle_mass")
        for field in fields:
            fv_hd = utils.generate_hash(
                self.field_values_test(ds_g8, field, particle_type=True)
            )
            hashes['field_values'][field] = fv_hd
        hashes = {'fields_g8' : hashes}
        utils.handle_hashes(self.save_dir, answer_file, hashes, self.answer_store)

    #-----
    # test_fields_g1
    #-----
    @utils.requires_ds(g1)
    def test_fields_g1(self, ds_g1):
        hashes = OrderedDict()
        hashes['field_values'] = OrderedDict()
        fields = ("particle_position_x", "particle_position_y",
                   "particle_position_z", "particle_mass")
        for field in fields:
            fv_hd = utils.generate_hash(
                self.field_values_test(ds_g1, field, particle_type=True)
            )
            hashes['field_values'][field] = fv_hd
        hashes = {'fields_g1' : hashes}
        utils.handle_hashes(self.save_dir, answer_file, hashes, self.answer_store)

    #-----
    # test_OWLSSubfindDataset
    #-----
    @requires_file(g1)
    def test_OWLSSubfindDataset(self, ds_g1):
        assert isinstance(ds_g1, OWLSSubfindDataset)
