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

from yt.testing import \
    assert_equal, \
    requires_file, \
    units_override_check
from yt.utilities.answer_testing.framework import \
    requires_ds, \
    small_patch_amr, \
    data_dir_load
from ..data_structures import FITSDataset

_fields_grs = ("temperature",)

grs = "radio_fits/grs-50-cube.fits"
@requires_ds(grs)
def test_grs():
    ds = data_dir_load(grs, cls=FITSDataset, kwargs={"nan_mask":0.0})
    yield assert_equal, str(ds), "grs-50-cube.fits"
    for test in small_patch_amr(ds, _fields_grs, input_center="c", input_weight="ones"):
        test_grs.__name__ = test.description
        yield test

_fields_vels = ("velocity_x","velocity_y","velocity_z")

vf = "UnigridData/velocity_field_20.fits"
@requires_ds(vf)
def test_velocity_field():
    ds = data_dir_load(vf, cls=FITSDataset)
    yield assert_equal, str(ds), "velocity_field_20.fits"
    for test in small_patch_amr(ds, _fields_vels, input_center="c", input_weight="ones"):
        test_velocity_field.__name__ = test.description
        yield test

@requires_file(vf)
def test_units_override():
    for test in units_override_check(vf):
        yield test

@requires_file(grs)
def test_FITSDataset():
    assert isinstance(data_dir_load(grs), FITSDataset)
