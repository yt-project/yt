"""
GDF frontend tests



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2016, yt Development Team.
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
from yt.frontends.gdf.api import GDFDataset

_fields = ("density", "velocity_x")

sedov = "sedov/sedov_tst_0004.h5"

@requires_ds(sedov)
def test_sedov_tunnel():
    ds = data_dir_load(sedov)
    yield assert_equal, str(ds), "sedov_tst_0004"
    for test in small_patch_amr(ds, _fields):
        test_sedov_tunnel.__name__ = test.description
        yield test


@requires_file(sedov)
def test_GDFDataset():
    assert isinstance(data_dir_load(sedov), GDFDataset)


@requires_file(sedov)
def test_units_override():
    for test in units_override_check(sedov):
        yield test
