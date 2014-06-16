"""
FITS frontend tests



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
    data_dir_load

_fields = ("intensity")

m33 = "radio_fits/m33_hi.fits"
@requires_ds(m33, big_data=True)
def test_m33():
    ds = data_dir_load(m33, nan_mask=0.0)
    yield assert_equal, str(ds), "m33_hi.fits"
    for test in small_patch_amr(m33, _fields):
        test_m33.__name__ = test.description
        yield test

_fields = ("temperature")

grs = "radio_fits/grs-50-cube.fits"
@requires_ds(grs)
def test_grs():
    ds = data_dir_load(grs, nan_mask=0.0)
    yield assert_equal, str(ds), "grs-50-cube.fits"
    for test in small_patch_amr(grs, _fields):
        test_grs.__name__ = test.description
        yield test

_fields = ("x-velocity","y-velocity","z-velocity")

vf = "UniformGrid/velocity_field_20.fits"
@requires_ds(vf)
def test_velocity_field():
    ds = data_dir_load(bf)
    yield assert_equal, str(ds), "velocity_field_20.fits"
    for test in small_patch_amr(vf, _fields):
        test_velocity_field.__name__ = test.description
        yield test
