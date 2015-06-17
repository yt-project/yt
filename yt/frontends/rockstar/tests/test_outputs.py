"""
Rockstar frontend tests using rockstar_halos dataset



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2015, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import os.path
from yt.testing import \
    assert_equal, \
    requires_file
from yt.utilities.answer_testing.framework import \
    FieldValuesTest, \
    requires_ds, \
    data_dir_load
from yt.frontends.rockstar.api import RockstarDataset

_fields = ("particle_position_x", "particle_position_y",
           "particle_position_z", "particle_mass")

r1 = "rockstar_halos/halos_0.0.bin"

@requires_ds(r1)
def test_fields_r1():
    ds = data_dir_load(r1)
    yield assert_equal, str(ds), os.path.basename(r1)
    for field in _fields:
        yield FieldValuesTest(r1, field, particle_type=True)

@requires_file(r1)
def test_RockstarDataset():
    assert isinstance(data_dir_load(r1), RockstarDataset)
