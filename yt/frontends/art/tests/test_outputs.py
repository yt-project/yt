"""
ART frontend tests using SFG1 a=0.330




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
from yt.frontends.art.api import ARTStaticOutput

_fields = ("Density", "particle_mass", ("all", "particle_position_x"))

sfg1 = "10MpcBox_csf512_a0.330.d"


@requires_pf(sfg1, big_data=True)
def test_sfg1():
    pf = data_dir_load(sfg1)
    yield assert_equal, str(pf), "10MpcBox_csf512_a0.330.d"
    dso = [None, ("sphere", ("max", (0.1, 'unitary')))]
    for field in _fields:
        for axis in [0, 1, 2]:
            for ds in dso:
                for weight_field in [None, "Density"]:
                    yield PixelizedProjectionValuesTest(
                        sfg1, axis, field, weight_field,
                        ds)
