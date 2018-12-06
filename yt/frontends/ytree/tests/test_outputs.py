"""
ytree frontend tests



"""

#-----------------------------------------------------------------------------
# Copyright (c) yt Development Team. All rights reserved.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import os.path
from yt.testing import \
    assert_equal, \
    requires_file, \
    requires_module
from yt.utilities.answer_testing.framework import \
    FieldValuesTest, \
    requires_ds, \
    data_dir_load
from yt.frontends.ytree.api import YTreeDataset

_fields = ("particle_position_x", "particle_position_y",
           "particle_position_z", "particle_mass")

arb = "arbor/arbor.h5"

@requires_ds(arb)
@requires_module('h5py')
def test_fields_arb():
    ds = data_dir_load(arb)
    assert_equal(str(ds), os.path.basename(arb))
    for field in _fields:
        yield FieldValuesTest(arb, field, particle_type=True)

@requires_file(arb)
@requires_module('h5py')
def test_RockstarDataset():
    assert isinstance(data_dir_load(arb), YTreeDataset)
