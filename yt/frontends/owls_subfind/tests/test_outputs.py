"""
Title: test_owls_subfind.py
Purpose: OWLSSubfind frontend tests using owls_fof_halos datasets
Notes:
    Copyright (c) 2013, yt Development Team.
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.
"""
from collections import OrderedDict

import pytest

from yt.frontends.owls_subfind.api import OWLSSubfindDataset
from yt.testing import requires_file
import yt.utilities.answer_testing.framework as fw
from yt.utilities.answer_testing import utils

# Test data
g1 = "owls_fof_halos/groups_001/group_001.0.hdf5"
g8 = "owls_fof_halos/groups_008/group_008.0.hdf5"


#============================================
#              TestOwlsSubfind
#============================================
@pytest.mark.skipif(not pytest.config.getvalue('--with-answer-testing'),
    reason="--with-answer-testing not set.")
@pytest.mark.usefixtures('answer_file')
class TestOwlsSubfind(fw.AnswerTest):
    #-----
    # test_fields_g8
    #-----
    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(g8)
    def test_fields_g8(self, ds_g8):
        self.hashes['field_values'] = OrderedDict()
        fields = ("particle_position_x", "particle_position_y",
                   "particle_position_z", "particle_mass")
        for field in fields:
            fv_hd = self.field_values_test(ds_g8, field, particle_type=True)
            self.hashes['field_values'][field] = fv_hd

    #-----
    # test_fields_g1
    #-----
    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(g1)
    def test_fields_g1(self, ds_g1):
        self.hashes['field_values'] = OrderedDict()
        fields = ("particle_position_x", "particle_position_y",
                   "particle_position_z", "particle_mass")
        for field in fields:
            fv_hd = self.field_values_test(ds_g1, field, particle_type=True)
            self.hashes['field_values'][field] = fv_hd

    #-----
    # test_OWLSSubfindDataset
    #-----
    @requires_file(g1)
    def test_OWLSSubfindDataset(self, ds_g1):
        assert isinstance(ds_g1, OWLSSubfindDataset)
