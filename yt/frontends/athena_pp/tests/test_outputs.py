"""
Athena++ frontend tests



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
from yt.testing import \
    assert_equal, \
    requires_file, \
    units_override_check, \
    assert_allclose
from yt.utilities.answer_testing.framework import \
    requires_ds, \
    small_patch_amr, \
    data_dir_load, \
    GenericArrayTest
from yt.frontends.athena_pp.api import AthenaPPDataset
from yt.convenience import load

_fields_disk = ("density", "velocity_r")

disk = "KeplerianDisk/disk.out1.00000.athdf"
@requires_ds(disk)
def test_disk():
    ds = data_dir_load(disk)
    assert_equal(str(ds), "disk.out1.00000")
    dd = ds.all_data()
    vol = (ds.domain_right_edge[0]**3-ds.domain_left_edge[0]**3)/3.0
    vol *= np.cos(ds.domain_left_edge[1])-np.cos(ds.domain_right_edge[1])
    vol *= ds.domain_right_edge[2].v-ds.domain_left_edge[2].v
    assert_allclose(dd.quantities.total_quantity("cell_volume"), vol)
    for field in _fields_disk:
        def field_func(name):
            return dd[field]
        yield GenericArrayTest(ds, field_func, args=[field])

_fields_AM06 = ("temperature", "density", "velocity_magnitude", "magnetic_field_x")

AM06 = "AM06/AM06.out1.00400.athdf"
@requires_ds(AM06)
def test_AM06():
    ds = data_dir_load(AM06)
    assert_equal(str(ds), "AM06.out1.00400")
    for test in small_patch_amr(ds, _fields_AM06):
        test_AM06.__name__ = test.description
        yield test

uo_AM06 = {
    'length_unit': (1.0, 'kpc'),
    'mass_unit': (1.0, 'Msun'),
    'time_unit': (1.0, 'Myr'),
}

@requires_file(AM06)
def test_AM06_override():
    # verify that overriding units causes derived unit values to be updated.
    # see issue #1259
    ds = load(AM06, units_override=uo_AM06)
    assert_equal(float(ds.magnetic_unit.in_units('gauss')),
                 9.01735778342523e-08)

@requires_file(AM06)
def test_units_override():
    units_override_check(AM06)

@requires_file(AM06)
def test_AthenaPPDataset():
    assert isinstance(data_dir_load(AM06), AthenaPPDataset)
