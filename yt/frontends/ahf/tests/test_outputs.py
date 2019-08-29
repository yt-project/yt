"""
Title: test_ahf.py
Purpose: Contains functions for answer testing the AHF frontend
Notes:
    Copyright (c) 2017, yt Development Team.
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.
"""
import os.path

from yt.frontends.ahf.api import AHFHalosDataset
from yt.testing import \
    assert_equal, \
    requires_file
from yt.utilities.answer_testing import \
    framework as fw, \
    utils

# Test data
ahf_halos = 'ahf_halos/snap_N64L16_135.parameter'


#============================================
#                   TestAHF
#============================================
class TestAHF(fw.AnswerTest):
    """
    Container for AHF frontend answer tests.

    Attributes:
    -----------
        pass

    Methods:
    --------
        pass
    """
    #-----
    # test_AHFHalosDataset
    #-----
    @requires_file(ahf_halos)
    def test_AHFHalosDataset(self, ds_ahf_halos):
        """
        Makes sure the dataset is loaded correctly.

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
        assert isinstance(ds_ahf_halos, AHFHalosDataset)

    #-----
    # test_fields_ahf_halos
    #-----
    @utils.requires_ds(ahf_halos)
    def test_fields_ahf_halos(self, ds_ahf_halos):
        """
        Runs the field_values_test on AHF fields.

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
        fields = ('particle_position_x', 'particle_position_y',
                   'particle_position_z', 'particle_mass')
        fv_hd = b''
        for field in fields:
            fv_hd += self.field_values_test(ds_ahf_halos, field, particle_type=True)
        # Hash the hex digest
        hashes = {'field_values' : utils.generate_hash(fv_hd)}
        # Save or compare answers
        utils.handle_hashes(self.save_dir, 'ahf', hashes, self.answer_store)
