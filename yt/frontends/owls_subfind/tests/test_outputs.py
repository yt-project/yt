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


#============================================
#              TestOwlsSubfind
#============================================
@pytest.mark.skipif(not pytest.config.getvalue('--with-answer-testing'),
    reason="--with-answer-testing not set.")
class TestOwlsSubfind(fw.AnswerTest):
    """
    Container for Owls subfind answer tests.

    Attributes:
    -----------
        pass

    Methods:
    --------
        pass
    """
    #-----
    # test_fields_g8
    #-----
    @utils.requires_ds(g8)
    def test_fields_g8(self, ds_g8):
        """
        Parameters:
        -----------
            pass

        Raises:
        -------
            pass

        Returns:
        --------
        """
        fv_hd = b''
        fields = ("particle_position_x", "particle_position_y",
                   "particle_position_z", "particle_mass")
        for field in fields:
            fv_hd += self.field_values_test(ds_g8, field, particle_type=True)
        hashes = {'field_values' : utils.generate_hash(fv_hd)}
        utils.handle_hashes(self.save_dir, 'owls-subfind-g8fields', hashes, self.answer_store)

    #-----
    # test_fields_g1
    #-----
    @utils.requires_ds(g1)
    def test_fields_g1(self, ds_g1):
        """
        Parameters:
        -----------
            pass

        Raises:
        -------
            pass

        Returns:
        --------
        """
        fv_hd = b''
        fields = ("particle_position_x", "particle_position_y",
                   "particle_position_z", "particle_mass")
        for field in fields:
            fv_hd += self.field_values_test(ds_g1, field, particle_type=True)
        hashes = {'field_values' : utils.generate_hash(fv_hd)}
        utils.handle_hashes(self.save_dir, 'owls-subfind-g1fields', hashes, self.answer_store)

    #-----
    # test_OWLSSubfindDataset
    #-----
    @requires_file(g1)
    def test_OWLSSubfindDataset(self, ds_g1):
        """
        Parameters:
        -----------
            pass

        Raises:
        -------
            pass

        Returns:
        --------
        """
        assert isinstance(ds_g1, OWLSSubfindDataset)
