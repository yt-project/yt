"""
GadgetFOF frontend tests using gadget_fof datasets



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
    requires_file, \
    assert_equal
from yt.utilities.answer_testing.framework import \
    FieldValuesTest, \
    requires_ds, \
    data_dir_load
from yt.frontends.gadget_fof.api import GadgetFOFDataset

_fields = ("particle_position_x", "particle_position_y",
           "particle_position_z", "particle_velocity_x",
           "particle_velocity_y", "particle_velocity_z",
           "particle_mass", "particle_identifier")

# a dataset with empty files
g5 = "gadget_fof_halos/groups_005/fof_subhalo_tab_005.0.hdf5"
g42 = "gadget_fof_halos/groups_042/fof_subhalo_tab_042.0.hdf5"


@requires_ds(g5)
def test_fields_g5():
    ds = data_dir_load(g5)
    yield assert_equal, str(ds), os.path.basename(g5)
    for field in _fields:
        yield FieldValuesTest(g5, field, particle_type=True)


@requires_ds(g42)
def test_fields_g42():
    ds = data_dir_load(g42)
    yield assert_equal, str(ds), os.path.basename(g42)
    for field in _fields:
        yield FieldValuesTest(g42, field, particle_type=True)

@requires_file(g42)
def test_GadgetFOFDataset():
    assert isinstance(data_dir_load(g42), GadgetFOFDataset)
