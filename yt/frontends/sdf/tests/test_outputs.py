"""
Title: test_sdf.py
Purpose: SDF frontend tests
Notes:
    Copyright (c) 2013, yt Development Team.
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.
"""
import socket

import numpy as np

from yt.extern.six.moves import urllib
from yt.frontends.sdf.api import SDFDataset
from yt.testing import \
    assert_equal, \
    requires_module
from yt.visualization.api import ProjectionPlot

import framework as fw
import utils

_fields = (('deposit', 'all_cic'))
scivis_data = "http://darksky.slac.stanford.edu/scivis2015/data/ds14_scivis_0128/ds14_scivis_0128_e4_dt04_1.0000"


#============================================
#                internet_on
#============================================
def internet_on():
    """
    Answer on http://stackoverflow.com/questions/3764291/checking-network-connection
    Better answer on:
    http://stackoverflow.com/questions/2712524/handling-urllib2s-timeout-python

    Parameters:
    -----------
        pass

    Raises:
    -------
        pass

    Returns:
    --------
        pass
    """
    try:
        urllib.request.urlopen(scivis_data, timeout=1)
        return True
    except urllib.error.URLError:
        return False
    except socket.timeout:
        return False


#============================================
#                  TestSDF
#============================================
class TestSDF(fw.AnswerTest):
    """
    Container for SDF frontend answer tests.

    Attributes:
    -----------
        pass

    Methods:
    --------
        pass
    """
    #-----
    # test_scivis
    #-----
    @requires_module('thingking')
    def test_scivis(self):
        if not internet_on():
            return
        ds = SDFDataset(scivis_data)
        assert_equal(str(ds), "ds14_scivis_0128_e4_dt04_1.0000")
        ad = ds.all_data()
        assert np.unique(ad['particle_position_x']).size > 1
        ProjectionPlot(ds, "z", _fields)
