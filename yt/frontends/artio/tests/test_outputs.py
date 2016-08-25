"""
ARTIO frontend tests




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.testing import \
    assert_equal, \
    requires_file, \
    units_override_check
from yt.utilities.answer_testing.framework import \
    requires_ds, \
    data_dir_load, \
    PixelizedProjectionValuesTest, \
    FieldValuesTest, \
    create_obj
from yt.frontends.artio.api import ARTIODataset

_fields = ("temperature", "density", "velocity_magnitude",
           ("deposit", "all_density"), ("deposit", "all_count"))

sizmbhloz = "sizmbhloz-clref04SNth-rs9_a0.9011/sizmbhloz-clref04SNth-rs9_a0.9011.art"
@requires_ds(sizmbhloz)
def test_sizmbhloz():
    ds = data_dir_load(sizmbhloz)
    ds.max_range = 1024*1024
    yield assert_equal, str(ds), "sizmbhloz-clref04SNth-rs9_a0.9011.art"
    dso = [ None, ("sphere", ("max", (0.1, 'unitary')))]
    for dobj_name in dso:
        for field in _fields:
            for axis in [0, 1, 2]:
                for weight_field in [None, "density"]:
                    yield PixelizedProjectionValuesTest(
                        ds, axis, field, weight_field,
                        dobj_name)
            yield FieldValuesTest(ds, field, dobj_name)
        dobj = create_obj(ds, dobj_name)
        s1 = dobj["ones"].sum()
        s2 = sum(mask.sum() for block, mask in dobj.blocks)
        yield assert_equal, s1, s2
    assert_equal(ds.particle_type_counts, {'N-BODY': 100000, 'STAR': 110650})


@requires_file(sizmbhloz)
def test_ARTIODataset():
    assert isinstance(data_dir_load(sizmbhloz), ARTIODataset)

@requires_file(sizmbhloz)
def test_units_override():
    for test in units_override_check(sizmbhloz):
        yield test
