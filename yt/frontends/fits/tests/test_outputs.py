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
from ..data_structures import FITSDataset

_fields_m33 = ("intensity",)

m33 = "radio_fits/m33_hi.fits"
@requires_pf(m33, big_data=True)
def test_m33():
    pf = data_dir_load(m33, cls=FITSDataset, kwargs={"nan_mask":0.0})
    yield assert_equal, str(pf), "m33_hi.fits"
    for test in small_patch_amr(m33, _fields_m33, input_center="c", input_weight="ones"):
        test_m33.__name__ = test.description
        yield test

_fields_grs = ("temperature",)

grs = "radio_fits/grs-50-cube.fits"
@requires_pf(grs)
def test_grs():
    pf = data_dir_load(grs, cls=FITSDataset, kwargs={"nan_mask":0.0})
    yield assert_equal, str(pf), "grs-50-cube.fits"
    for test in small_patch_amr(grs, _fields_grs, input_center="c", input_weight="ones"):
        test_grs.__name__ = test.description
        yield test

_fields_vels = ("x-velocity","y-velocity","z-velocity")

vf = "UniformGrid/velocity_field_20.fits"
@requires_pf(vf)
def test_velocity_field():
    pf = data_dir_load(bf, cls=FITSDataset)
    yield assert_equal, str(pf), "velocity_field_20.fits"
    for test in small_patch_amr(vf, _fields_vels):
        test_velocity_field.__name__ = test.description
        yield test
