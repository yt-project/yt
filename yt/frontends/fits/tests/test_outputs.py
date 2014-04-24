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
    requires_pf, \
    small_patch_amr, \
    big_patch_amr, \
    data_dir_load

_fields = ("intensity")

m33 = "fits/m33_hi.fits"
@requires_pf(m33, big_data=True)
def test_m33():
    pf = data_dir_load(m33)
    yield assert_equal, str(pf), "m33_hi.fits"
    for test in small_patch_amr(m33, _fields):
        test_m33.__name__ = test.description
        yield test

_fields = ("x-velocity","y-velocity","z-velocity")

vf = "fits/velocity_field_20.fits"
@requires_pf(vf)
def test_velocity_field():
    pf = data_dir_load(bf)
    yield assert_equal, str(pf), "velocity_field_20.fits"
    for test in small_patch_amr(vf, _fields):
        test_velocity_field.__name__ = test.description
        yield test
