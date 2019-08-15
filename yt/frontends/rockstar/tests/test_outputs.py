"""
Title: test_rockstar.py
Purpose: Rockstar frontend tests using rockstar_halos dataset
Notes:
    Copyright (c) 2015, yt Development Team.
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.
"""
import os.path

from yt.testing import \
    assert_equal, \
    requires_file
from yt.frontends.rockstar.api import RockstarDataset

import framework as fw
import utils


# Test data
r1 = "rockstar_halos/halos_0.0.bin"


# Globals
_fields = ("particle_position_x", "particle_position_y",
           "particle_position_z", "particle_mass")


#============================================
#               TestRockstar
#============================================
class TestRockstar(fw.AnswerTest):
    """
    Container for rockstar frontend answer tests.

    Attributes:
    -----------
        pass

    Methods:
    --------
        pass
    """
    #-----
    # test_fields_r1
    #-----
    @utils.requires_ds(r1)
    def test_fields_r1(self):
        """
        Parameters:
        -----------
            pass

        Raises:
        -------
            pass

        Returns:
        --------
            pass
        """
        fv_hd = b''
        ds = utils.data_dir_load(r1)
        assert_equal(str(ds), os.path.basename(r1))
        for field in _fields:
            fv_hd += self.field_values_test(ds, field, particle_type=True)
        hashes = {'field_values' : utils.generate_hash(fv_hd)}
        utils.handle_hashes(self.save_dir, 'rockstar-test-fields-r1', hashes, self.answer_store)

    #-----
    # test_RockstarDataset
    #-----
    @requires_file(r1)
    def test_RockstarDataset(self):
        """
        Parameters:
        -----------
            pass

        Raises:
        -------
            pass

        Returns:
        --------
            pass
        """
        assert isinstance(utils.data_dir_load(r1), RockstarDataset)
