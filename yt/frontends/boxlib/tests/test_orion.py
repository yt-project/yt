"""
Orion frontend tests



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
from yt.frontends.boxlib.api import OrionDataset

# We don't do anything needing ghost zone generation right now, because these
# are non-periodic datasets.
_fields = ("temperature", "density", "velocity_magnitude")

radadvect = "RadAdvect/plt00000"
@requires_ds(radadvect)
def test_radadvect():
    ds = data_dir_load(radadvect)
    yield assert_equal, str(ds), "plt00000"
    for test in small_patch_amr(ds, _fields):
        test_radadvect.__name__ = test.description
        yield test

rt = "RadTube/plt00500"
@requires_ds(rt)
def test_radtube():
    ds = data_dir_load(rt)
    yield assert_equal, str(ds), "plt00500"
    for test in small_patch_amr(ds, _fields):
        test_radtube.__name__ = test.description
        yield test

star = "StarParticles/plrd01000"
@requires_ds(star)
def test_star():
    ds = data_dir_load(star)
    yield assert_equal, str(ds), "plrd01000"
    for test in small_patch_amr(ds, _fields):
        test_star.__name__ = test.description
        yield test

@requires_file(rt)
def test_OrionDataset():
    assert isinstance(data_dir_load(rt), OrionDataset)

@requires_file(rt)
def test_units_override():
    for test in units_override_check(rt):
        yield test

