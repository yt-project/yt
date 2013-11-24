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

from yt.testing import *
from yt.utilities.answer_testing.framework import \
    requires_pf, \
    data_dir_load, \
    PixelizedProjectionValuesTest, \
    FieldValuesTest, \
    create_obj
from yt.frontends.artio.api import ARTIOStaticOutput

_fields = ("Temperature", "Density", "VelocityMagnitude",
           ("deposit", "all_density"), ("deposit", "all_count")) 

sizmbhloz = "sizmbhloz-clref04SNth-rs9_a0.9011/sizmbhloz-clref04SNth-rs9_a0.9011.art"
@requires_pf(sizmbhloz)
def test_sizmbhloz():
    pf = data_dir_load(sizmbhloz)
    pf.max_range = 1024*1024
    yield assert_equal, str(pf), "sizmbhloz-clref04SNth-rs9_a0.9011.art"
    dso = [ None, ("sphere", ("max", (0.1, 'unitary')))]
    for ds in dso:
        for field in _fields:
            for axis in [0, 1, 2]:
                for weight_field in [None, "Density"]:
                    yield PixelizedProjectionValuesTest(
                        sizmbhloz, axis, field, weight_field,
                        ds)
            yield FieldValuesTest(sizmbhloz, field, ds)
        dobj = create_obj(pf, ds)
        s1 = dobj["Ones"].sum()
        s2 = sum(mask.sum() for block, mask in dobj.blocks)
        yield assert_equal, s1, s2
