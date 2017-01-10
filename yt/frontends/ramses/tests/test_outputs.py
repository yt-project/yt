"""
RAMSES frontend tests




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
from yt.frontends.ramses.api import RAMSESDataset
import os
import yt

_fields = ("temperature", "density", "velocity_magnitude",
           ("deposit", "all_density"), ("deposit", "all_count"))

output_00080 = "output_00080/info_00080.txt"
@requires_ds(output_00080)
def test_output_00080():
    ds = data_dir_load(output_00080)
    yield assert_equal, str(ds), "info_00080"
    dso = [ None, ("sphere", ("max", (0.1, 'unitary')))]
    for dobj_name in dso:
        for field in _fields:
            for axis in [0, 1, 2]:
                for weight_field in [None, "density"]:
                    yield PixelizedProjectionValuesTest(
                        output_00080, axis, field, weight_field,
                        dobj_name)
            yield FieldValuesTest(output_00080, field, dobj_name)
        dobj = create_obj(ds, dobj_name)
        s1 = dobj["ones"].sum()
        s2 = sum(mask.sum() for block, mask in dobj.blocks)
        yield assert_equal, s1, s2
    assert_equal(ds.particle_type_counts, {'io': 1090895})

@requires_file(output_00080)
def test_RAMSESDataset():
    assert isinstance(data_dir_load(output_00080), RAMSESDataset)

@requires_file(output_00080)
def test_units_override():
    for test in units_override_check(output_00080):
        yield test


ramsesNonCosmo = 'DICEGalaxyDisk_nonCosmological/output_00002'
@requires_file(ramsesNonCosmo)
def test_unit_non_cosmo():
    ds = yt.load(os.path.join(ramsesNonCosmo, 'info_00002.txt'))

    expected_raw_time = 0.0299468077820411 # in ramses unit
    yield assert_equal, ds.current_time.value, expected_raw_time

    expected_time = 14087886140997.336 # in seconds
    assert_equal(ds.current_time.in_units('s').value, expected_time)
