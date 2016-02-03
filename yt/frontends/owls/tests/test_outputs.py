"""
OWLS frontend tests using the snapshot_033 dataset




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from collections import OrderedDict

from yt.testing import \
    requires_file
from yt.utilities.answer_testing.framework import \
    requires_ds, \
    data_dir_load, \
    sph_answer
from yt.frontends.owls.api import OWLSDataset

os33 = "snapshot_033/snap_033.0.hdf5"

# This maps from field names to weight field names to use for projections
_fields = OrderedDict(
    [
        (("gas", "density"), None),
        (("gas", "temperature"), None),
        (("gas", "temperature"), ("gas", "density")),
        (('gas', 'He_p0_number_density'), None),
        (('gas', 'velocity_magnitude'), None),
        (("deposit", "all_density"), None),
        (("deposit", "all_count"), None),
        (("deposit", "all_cic"), None),
        (("deposit", "PartType0_density"), None),
        (("deposit", "PartType4_density"), None),
    ]
)


@requires_ds(os33, big_data=True)
def test_snapshot_033():
    ds = data_dir_load(os33)
    for test in sph_answer(ds, 'snap_033', 2*128**3, _fields):
        test_snapshot_033.__name__ = test.description
        yield test


@requires_file(os33)
def test_OWLSDataset():
    assert isinstance(data_dir_load(os33), OWLSDataset)
