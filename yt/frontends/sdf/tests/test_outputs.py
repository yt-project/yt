"""
SDF frontend tests

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013-2020, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
import socket
from yt.testing import assert_equal
from yt.frontends.sdf.api import SDFDataset
from yt.visualization.api import ProjectionPlot
from yt.testing import requires_module
from yt.extern.six.moves import urllib


_fields = (('deposit', 'all_cic'))
slac_scivis_data = "http://darksky.slac.stanford.edu/scivis2015/data/ds14_scivis_0128/ds14_scivis_0128_e4_dt04_1.0000"
ncsa_scivis_data = "http://use.yt/upload/744abba3"
scivis_data = ncsa_scivis_data

# Answer on http://stackoverflow.com/questions/3764291/checking-network-connection
# Better answer on http://stackoverflow.com/questions/2712524/handling-urllib2s-timeout-python
def internet_on():
    try:
        urllib.request.urlopen(scivis_data, timeout=1)
        return True
    except urllib.error.URLError:
        return False
    except socket.timeout:
        return False


@requires_module('thingking')
def test_scivis():
    if not internet_on():
        return
    return # HOTFIX: See discussion in 2334
    ds = SDFDataset(scivis_data)
    if scivis_data == slac_scivis_data:
        assert_equal(str(ds), "ds14_scivis_0128_e4_dt04_1.0000")
    else:
        assert_equal(str(ds), "744abba3")
    ad = ds.all_data()
    assert np.unique(ad['particle_position_x']).size > 1
    ProjectionPlot(ds, "z", _fields)
