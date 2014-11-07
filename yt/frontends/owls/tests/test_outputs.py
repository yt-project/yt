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
    requires_ds, \
    small_patch_amr, \
    big_patch_amr, \
    data_dir_load, \
    PixelizedProjectionValuesTest, \
    FieldValuesTest, \
    create_obj
from yt.frontends.owls.api import OWLSDataset

_fields = (("deposit", "all_density"), ("deposit", "all_count"),
           ("deposit", "PartType0_density"),
           ("deposit", "PartType4_density"))

os33 = "snapshot_033/snap_033.0.hdf5"
@requires_ds(os33, big_data=True)
def test_snapshot_033():
    ds = data_dir_load(os33)
    yield assert_equal, str(ds), "snap_033"
    dso = [ None, ("sphere", ("c", (0.1, 'unitary')))]
    dd = ds.all_data()
    yield assert_equal, dd["particle_position"].shape[0], 2*(128*128*128)
    yield assert_equal, dd["particle_position"].shape[1], 3
    tot = sum(dd[ptype,"particle_position"].shape[0]
              for ptype in ds.particle_types if ptype != "all")
    yield assert_equal, tot, (2*128*128*128)
    for dobj_name in dso:
        for field in _fields:
            for axis in [0, 1, 2]:
                for weight_field in [None, "density"]:
                    yield PixelizedProjectionValuesTest(
                        os33, axis, field, weight_field,
                        dobj_name)
            yield FieldValuesTest(os33, field, dobj_name)
        dobj = create_obj(ds, dobj_name)
        s1 = dobj["ones"].sum()
        s2 = sum(mask.sum() for block, mask in dobj.blocks)
        yield assert_equal, s1, s2


@requires_file(os33)
def test_OWLSDataset():
    assert isinstance(data_dir_load(os33), OWLSDataset)
