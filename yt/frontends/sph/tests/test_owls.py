"""
Tipsy tests using the OWLS HDF5-Gadget dataset




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
    requires_pf, \
    small_patch_amr, \
    big_patch_amr, \
    data_dir_load, \
    PixelizedProjectionValuesTest, \
    FieldValuesTest
from yt.frontends.sph.api import OWLSDataset

_fields = (("deposit", "all_density"), ("deposit", "all_count"),
           ("deposit", "PartType0_density"),
           ("deposit", "PartType4_density"))

os33 = "snapshot_033/snap_033.0.hdf5"
@requires_pf(os33)
def test_snapshot_033():
    pf = data_dir_load(os33)
    yield assert_equal, str(pf), "snap_033"
    dso = [ None, ("sphere", ("c", (0.1, 'unitary')))]
    dd = pf.h.all_data()
    yield assert_equal, dd["Coordinates"].shape[0], 2*(128*128*128)
    yield assert_equal, dd["Coordinates"].shape[1], 3
    tot = sum(dd[ptype,"Coordinates"].shape[0]
              for ptype in pf.particle_types if ptype != "all")
    yield assert_equal, tot, (2*128*128*128)
    for field in _fields:
        for axis in [0, 1, 2]:
            for ds in dso:
                for weight_field in [None, "Density"]:
                    yield PixelizedProjectionValuesTest(
                        os33, axis, field, weight_field,
                        ds)
                yield FieldValuesTest(
                        os33, field, ds)

