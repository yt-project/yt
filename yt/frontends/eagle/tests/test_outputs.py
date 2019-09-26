"""
Title: test_eagle.py
Purpose: Eagle frontend tests using the snapshot_028_z000p000 dataset
    Copyright (c) 2015, yt Development Team.
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.
"""
import pytest

from yt.testing import requires_file
from yt.frontends.eagle.api import EagleDataset
import yt.utilities.answer_testing.framework as fw
from yt.utilities.answer_testing import utils


# Test data
s28 = "snapshot_028_z000p000/snap_028_z000p000.0.hdf5"
s399 = "snipshot_399_z000p000/snip_399_z000p000.0.hdf5"


#============================================
#                 TestEagle
#============================================
@pytest.mark.skipif(not pytest.config.getvalue('--with-answer-testing'),
    reason="--with-answer-testing not set.")
@pytest.mark.usefixtures('answer_file')
class TestEagle(fw.AnswerTest):
    #-----
    # test_EagleDataset
    #-----
    @requires_file(s28)
    def test_EagleDataset(self):
        assert isinstance(utils.data_dir_load(s28), EagleDataset)

    #-----
    # test_Snipshot
    #-----
    @requires_file(s399)
    def test_Snipshot(self):
        ds = utils.data_dir_load(s399)
        ds.index
        assert isinstance(ds, EagleDataset)
