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

from yt.testing import \
    requires_file
from yt.utilities.answer_testing.framework import \
    requires_ds, \
    data_dir_load, \
    sph_answer_test
from yt.frontends.owls.api import OWLSDataset

os33 = "snapshot_033/snap_033.0.hdf5"

_fields = (
    ("gas", "density")
    ("gas", "temperature"),
    ('gas', 'He_p0_number_density'),
    ('gas', 'N_p1_number_density'),
    ('gas', 'velocity_magnitude'),
    ("deposit", "all_density"),
    ("deposit", "all_count"),
    ("deposit", "all_cic")
    ("deposit", "PartType0_density"),
    ("deposit", "PartType4_density"))


@requires_ds(os33, big_data=True)
def test_snapshot_033():
    yield sph_answer_test(os33, 'snap_033', 2*128**3, _fields)


@requires_file(os33)
def test_OWLSDataset():
    assert isinstance(data_dir_load(os33), OWLSDataset)
