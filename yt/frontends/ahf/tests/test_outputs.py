"""
AHF frontend tests using ahf_halos dataset



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2017, yt Development Team.
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
from yt.frontends.ahf.api import AHFHalosDataset

_fields = ('particle_position_x', 'particle_position_y',
           'particle_position_z', 'particle_mass')

ahf_halos = 'ahf_halos/snap_N64L16_135.parameter'


@requires_ds(ahf_halos)
def test_fields_sample():
    ds = data_dir_load(ahf_halos)
    assert_equal(str(ahf_halos), os.path.basename(ahf_halos))
    for field in _fields:
        yield FieldValuesTest(ahf_halos, field, particle_type=True)

@requires_file(ahf_halos)
def test_AHFHalosDataset():
    assert isinstance(data_dir_load(ahf_halos), AHFHalosDataset)
