"""
SDF frontend tests

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
from yt.visualization.api import ProjectionPlot

_fields = (('deposit','all_cic'))

import urllib2

scivis_data = "http://darksky.slac.stanford.edu/scivis2015/data/ds14_scivis_0128/ds14_scivis_0128_e4_dt04_1.0000"

# Answer on http://stackoverflow.com/questions/3764291/checking-network-connection
def internet_on():
    try:
        response=urllib2.urlopen(scivis_data,timeout=1)
        return True
    except urllib2.URLError as err: pass
    return False

def test_scivis():
    if not internet_on(): return
    ds = SDFDataset(scivis_data)
    yield assert_equal, str(ds), "ds14_scivis_0128_e4_dt04_1.0000"
    ad = ds.all_data()
    assert np.unique(ad['particle_position_x']).size > 1
    p = ProjectionPlot(ds, "z", _fields)
