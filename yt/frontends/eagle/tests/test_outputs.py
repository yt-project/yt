"""
Eagle frontend tests using the snapshot_028_z000p000 dataset




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2015, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.testing import \
    requires_file
from yt.utilities.answer_testing.framework import \
    data_dir_load
from yt.frontends.eagle.api import EagleDataset

s28 = "snapshot_028_z000p000/snap_028_z000p000.0.hdf5"
@requires_file(s28)
def test_EagleDataset():
    assert isinstance(data_dir_load(s28), EagleDataset)

s399 = "snipshot_399_z000p000/snip_399_z000p000.0.hdf5"
@requires_file(s399)
def test_Snipshot():
    ds = data_dir_load(s399)
    ds.index
    assert isinstance(ds, EagleDataset)
