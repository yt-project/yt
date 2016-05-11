"""
Gizmo frontend tests




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2015, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from collections import OrderedDict

from yt.utilities.answer_testing.framework import \
    data_dir_load, \
    requires_ds, \
    sph_answer
from yt.frontends.gizmo.api import GizmoDataset

FIRE_m12i = 'FIRE_M12i_ref11/snapshot_600.hdf5'

# This maps from field names to weight field names to use for projections
fields = OrderedDict(
    [
        (("gas", "density"), None),
        (("gas", "temperature"), None),
        (("gas", "temperature"), ('gas', 'density')),
        (('gas', 'velocity_magnitude'), None),
        (("deposit", "all_density"), None),
        (("deposit", "all_count"), None),
        (("deposit", "all_cic"), None),
        (("deposit", "PartType0_density"), None),
    ]
)

@requires_ds(FIRE_m12i)
def test_GizmoDataset():
    ds = data_dir_load(FIRE_m12i)
    assert isinstance(ds, FIREDataset)
    for test in sph_answer(ds, 'FIRE_m12i', 4786950, fields):
        test_GizmoDataset.__name__ = test.description
        yield test
    assert False
