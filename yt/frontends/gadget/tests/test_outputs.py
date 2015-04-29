"""
Gadget frontend tests using the IsothermalCollapse dataset




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2015, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.testing import *
from yt.utilities.answer_testing.framework import \
    requires_ds, \
    data_dir_load
from yt.frontends.gadget.api import GadgetHDF5Dataset

isothermal = "IsothermalCollapse/snap_505.hdf5"
@requires_file(isothermal)
def test_GadgetDataset():
    assert isinstance(data_dir_load(isothermal), GadgetHDF5Dataset)
