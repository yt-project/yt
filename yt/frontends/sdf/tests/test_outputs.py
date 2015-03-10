"""
Enzo frontend tests using moving7



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.testing import *
import numpy as np
from yt.frontends.sdf.api import SDFDataset

_fields = (('deposit','all_cic'))

def test_scivis():
    ds = SDFDataset("http://darksky.slac.stanford.edu/scivis2015/data/ds14_scivis_0128/ds14_scivis_0128_e4_dt04_1.0000")
    yield assert_equal, str(ds), "ds14_scivis_0128_e4_dt04_1.0000"
    ad = ds.all_data()
    assert np.unique(ad['particle_position_x']).size > 1
