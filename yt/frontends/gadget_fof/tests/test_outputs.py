"""
GadgetFOF frontend tests using owls_fof_halos datasets



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import os.path
from yt.testing import \
    assert_equal
from yt.utilities.answer_testing.framework import \
    FieldValuesTest, \
    requires_ds, \
    data_dir_load
from yt.frontends.gadget_fof.api import GadgetFOFDataset

_fields = ("particle_position_x", "particle_position_y",
           "particle_position_z", "particle_mass")

# a dataset with empty files
g1 = "" # TBD
g8 = "" # TBD


@requires_ds(g8)
def test_fields_g8():
    ds = data_dir_load(g8)
    yield assert_equal, str(ds), os.path.basename(g8)
    for field in _fields:
        yield FieldValuesTest(g8, field, particle_type=True)


@requires_ds(g1)
def test_fields_g1():
    ds = data_dir_load(g1)
    yield assert_equal, str(ds), os.path.basename(g1)
    for field in _fields:
        yield FieldValuesTest(g1, field, particle_type=True)

@requires_file(g1)
def test_GadgetFOFDataset():
    assert isinstance(data_dir_load(g1), GadgetFOFDataset)
