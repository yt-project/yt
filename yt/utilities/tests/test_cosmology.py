"""
Test cosmology calculator.




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2014, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from yt.testing import assert_rel_equal
from yt.utilities.cosmology import \
     Cosmology

def test_hubble_time():
    """
    Make sure hubble_time and t_from_z functions agree.

    """

    for i in range(10):
        co = Cosmology()
        # random sample over interval (-1,100]
        z = -101 * np.random.random() + 100
        yield assert_rel_equal, co.hubble_time(z), co.t_from_z(z), 5

def test_z_t_conversion():
    """
    Make sure t_from_z and z_from_t are consistent.

    """

    for i in range(10):
        co = Cosmology()
        # random sample over interval (-1,100]
        z1 = -101 * np.random.random() + 100
        t = co.t_from_z(z1)
        z2 = co.z_from_t(t)
        yield assert_rel_equal, z1, z2, 10
