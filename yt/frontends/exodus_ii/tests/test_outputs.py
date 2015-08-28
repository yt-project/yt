"""
Exodus II frontend tests



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
    data_dir_load
from yt.frontends.exodus_ii.api import ExodusIIDataset

_fields = ("temperature", "density", "velocity_magnitude",
           "velocity_divergence")

print "Is this thing on"

out00 = "ExodusII/out.e"
@requires_ds(out00)
def test_printout():
    print "Is it me you're looking for?"
    ds = data_dir_load(out00)
    yield assert_equal, str(ds), "out.e"


