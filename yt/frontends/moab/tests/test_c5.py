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

from yt.testing import *
from yt.utilities.answer_testing.framework import \
    requires_pf, \
    small_patch_amr, \
    big_patch_amr, \
    data_dir_load, \
    PixelizedProjectionValuesTest, \
    FieldValuesTest
from yt.frontends.moab.api import MoabHex8Dataset

_fields = (("gas", "flux"),
          )

c5 = "c5/c5.h5m"
@requires_pf(c5)
def test_cantor_5():
    np.random.seed(0x4d3d3d3)
    pf = data_dir_load(c5)
    yield assert_equal, str(pf), "c5"
    dso = [ None, ("sphere", ("c", (0.1, 'unitary'))),
                  ("sphere", ("c", (0.2, 'unitary')))]
    dd = pf.h.all_data()
    yield assert_almost_equal, pf.h.get_smallest_dx(), 0.00411522633744843, 10
    yield assert_equal, dd["x"].shape[0], 63*63*63
    yield assert_almost_equal, dd["CellVolumeCode"].sum(dtype="float64"), 1.0, 10
    for offset_1 in [1e-9, 1e-4, 0.1]:
        for offset_2 in [1e-9, 1e-4, 0.1]:
            ray = pf.h.ray(pf.domain_left_edge + offset_1,
                           pf.domain_right_edge - offset_2)
            yield assert_almost_equal, ray["dts"].sum(dtype="float64"), 1.0, 8
    for i, p1 in enumerate(np.random.random((5, 3))):
        for j, p2 in enumerate(np.random.random((5, 3))):
            ray = pf.h.ray(p1, p2)
            yield assert_almost_equal, ray["dts"].sum(dtype="float64"), 1.0, 8
    for field in _fields:
        for axis in [0, 1, 2]:
            for ds in dso:
                yield FieldValuesTest(c5, field, ds)

