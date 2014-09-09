"""
Chombo frontend tests



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
    assert_equal
from yt.utilities.answer_testing.framework import \
    requires_ds, \
    small_patch_amr, \
    data_dir_load
from yt.frontends.chombo.api import \
    ChomboDataset, \
    Orion2Dataset

_fields = ("density", "velocity_magnitude",  # "velocity_divergence",
           "magnetic_field_x")

gc = "GaussianCloud/data.0077.3d.hdf5"
@requires_ds(gc)
def test_gc():
    ds = data_dir_load(gc)
    yield assert_equal, str(ds), "data.0077.3d.hdf5"
    for test in small_patch_amr(gc, _fields):
        test_gc.__name__ = test.description
        yield test

tb = "TurbBoxLowRes/data.0005.3d.hdf5"
@requires_ds(tb)
def test_tb():
    ds = data_dir_load(tb)
    yield assert_equal, str(ds), "data.0005.3d.hdf5"
    for test in small_patch_amr(tb, _fields):
        test_tb.__name__ = test.description
        yield test

_zp_fields = ("rhs", "phi", "gravitational_field_x",
              "gravitational_field_y")
zp = "ZeldovichPancake/plt32.2d.hdf5"
@requires_ds(zp)
def test_zp():
    ds = data_dir_load(zp)
    yield assert_equal, str(ds), "plt32.2d.hdf5"
    for test in small_patch_amr(zp, _zp_fields, input_center="c",
                                input_weight="rhs"):
        test_tb.__name__ = test.description
        yield test


@requires_file(zp)
def test_ChomboDataset():
    assert isinstance(data_dir_load(zp), ChomboDataset)


@requires_file(gc)
def test_Orion2Dataset():
    assert isinstance(data_dir_load(gc), Orion2Dataset)


#@requires_file(kho)
#def test_PlutoDataset():
#    assert isinstance(data_dir_load(kho), PlutoDataset)
