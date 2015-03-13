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
from yt.testing import requires_module

_fields = (('deposit','all_cic'))

import urllib2
import socket

scivis_data = "http://darksky.slac.stanford.edu/scivis2015/data/ds14_scivis_0128/ds14_scivis_0128_e4_dt04_1.0000"

# Answer on http://stackoverflow.com/questions/3764291/checking-network-connection
# Better answer on http://stackoverflow.com/questions/2712524/handling-urllib2s-timeout-python
def internet_on():
    try:
        urllib2.urlopen("http://www.example.com", timeout = 1)
        return True
    except urllib2.URLError, e:
        # For Python 2.6
        if isinstance(e.reason, socket.timeout):
            return False
        else:
            # reraise the original error
            raise
    except socket.timeout, e:
        # For Python 2.7
        return False

@requires_module('thingking')
def test_scivis():
    if not internet_on(): return
    ds = SDFDataset(scivis_data)
    yield assert_equal, str(ds), "ds14_scivis_0128_e4_dt04_1.0000"
    ad = ds.all_data()
    assert np.unique(ad['particle_position_x']).size > 1
    p = ProjectionPlot(ds, "z", _fields)
