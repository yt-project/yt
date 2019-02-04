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

def test_z_t_conversion():
    """
    Make sure t_from_z and z_from_t are consistent.

    """

    co = Cosmology()
    # random sample in log(a) from -6 to 6
    my_random = np.random.RandomState(2419)
    la = 12 * my_random.random_sample(10000) - 6
    z1 = 1 / np.power(10, la) - 1
    t = co.t_from_z(z1)
    z2 = co.z_from_t(t)
    assert_rel_equal(z1, z2, 4)
