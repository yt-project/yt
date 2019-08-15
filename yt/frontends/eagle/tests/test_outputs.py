"""
Title: test_eagle.py
Purpose: Eagle frontend tests using the snapshot_028_z000p000 dataset
    Copyright (c) 2015, yt Development Team.
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.
"""
from yt.testing import requires_file
from yt.frontends.eagle.api import EagleDataset

import framework as fw
import utils


# Test data
s28 = "snapshot_028_z000p000/snap_028_z000p000.0.hdf5"
s399 = "snipshot_399_z000p000/snip_399_z000p000.0.hdf5"


#============================================
#                 TestEagle
#============================================
class TestEagle(fw.AnswerTest):
    """
    Container for eagle frontend answer tests.

    Attributes:
    -----------
        pass

    Methods:
    --------
        pass
    """
    #-----
    # test_EagleDataset
    #-----
    @requires_file(s28)
    def test_EagleDataset(self):
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
        assert isinstance(utils.data_dir_load(s28), EagleDataset)

    #-----
    # test_Snipshot
    #-----
    @requires_file(s399)
    def test_Snipshot(self):
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
        ds = utils.data_dir_load(s399)
        ds.index
        assert isinstance(ds, EagleDataset)
