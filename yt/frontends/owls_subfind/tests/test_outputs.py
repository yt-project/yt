"""
OWLSSubfind frontend tests using owls_fof_halos datasets



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.testing import *
from yt.utilities.answer_testing.framework import \
    FieldValuesTest, \
    requires_ds, \
    data_dir_load

_fields = ("particle_position_x", "particle_position_y",
           "particle_position_z", "particle_mass")

g8 = "owls_fof_halos/groups_008/group_008.0.hdf5"
@requires_ds(g8)
def test_fields_g8():
    ds = data_dir_load(g8)
    yield assert_equal, str(ds), "group_008.0.hdf5"
    for field in _fields:
        yield FieldValuesTest(g8, field, particle_type=True)

# a dataset with empty files
g3 = "owls_fof_halos/groups_003/group_003.0.hdf5"
@requires_ds(g3)
def test_fields_g3():
    ds = data_dir_load(g3)
    yield assert_equal, str(ds), "group_003.0.hdf5"
    for field in _fields:
        yield FieldValuesTest(g3, field, particle_type=True)
