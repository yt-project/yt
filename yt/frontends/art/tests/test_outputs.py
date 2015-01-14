"""
ART frontend tests using D9p a=0.500




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.testing import \
    requires_file, \
    assert_equal, \
    units_override_check
from yt.utilities.answer_testing.framework import \
    requires_ds, \
    big_patch_amr, \
    PixelizedProjectionValuesTest, \
    data_dir_load
from yt.frontends.art.api import ARTDataset

_fields = ("density", "temperature", "particle_mass", ("all", "particle_position_x"))

d9p = "D9p_500/10MpcBox_HartGal_csf_a0.500.d"

@requires_ds(d9p, big_data=True)
def test_d9p():
    ds = data_dir_load(d9p)
    yield assert_equal, str(ds), "10MpcBox_HartGal_csf_a0.500.d"
    for test in big_patch_amr(d9p, _fields):
        test_d9p.__name__ = test.description
        yield test
    dso = [None, ("sphere", ("max", (0.1, 'unitary')))]
    for field in _fields:
        for axis in [0, 1, 2]:
            for dobj_name in dso:
                for weight_field in [None, "density"]:
                    yield PixelizedProjectionValuesTest(
                        d9p, axis, field, weight_field,
                        dobj_name)


@requires_file(d9p)
def test_ARTDataset():
    assert isinstance(data_dir_load(d9p), ARTDataset)

@requires_file(d9p)
def test_units_override():
    for test in units_override_check(d9p):
        yield test

