"""
Tests of semi-structured meshes in MoabHex8 format.




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
    assert_almost_equal, \
    requires_file, \
    units_override_check
from yt.utilities.answer_testing.framework import \
    requires_ds, \
    data_dir_load, \
    FieldValuesTest
from yt.frontends.moab.api import MoabHex8Dataset

_fields = (("moab", "flux"),
          )

c5 = "c5/c5.h5m"
@requires_ds(c5)
def test_cantor_5():
    np.random.seed(0x4d3d3d3)
    ds = data_dir_load(c5)
    yield assert_equal, str(ds), "c5"
    dso = [ None, ("sphere", ("c", (0.1, 'unitary'))),
                  ("sphere", ("c", (0.2, 'unitary')))]
    dd = ds.all_data()
    yield assert_almost_equal, ds.index.get_smallest_dx(), 0.00411522633744843, 10
    yield assert_equal, dd["x"].shape[0], 63*63*63
    yield assert_almost_equal, \
        dd["cell_volume"].in_units("code_length**3").sum(dtype="float64").d, \
        1.0, 10
    for offset_1 in [1e-9, 1e-4, 0.1]:
        for offset_2 in [1e-9, 1e-4, 0.1]:
            DLE = ds.domain_left_edge
            DRE = ds.domain_right_edge
            ray = ds.ray(DLE + offset_1 * DLE.uq,
                         DRE - offset_2 * DRE.uq)
            yield assert_almost_equal, ray["dts"].sum(dtype="float64"), 1.0, 8
    for i, p1 in enumerate(np.random.random((5, 3))):
        for j, p2 in enumerate(np.random.random((5, 3))):
            ray = ds.ray(p1, p2)
            yield assert_almost_equal, ray["dts"].sum(dtype="float64"), 1.0, 8
    for field in _fields:
        for dobj_name in dso:
            yield FieldValuesTest(c5, field, dobj_name)


@requires_file(c5)
def test_MoabHex8Dataset():
    assert isinstance(data_dir_load(c5), MoabHex8Dataset)

@requires_file(c5)
def test_units_override():
    for test in units_override_check(c5):
        yield test
