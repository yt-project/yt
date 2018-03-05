"""Tests for selection_routines."""

#-----------------------------------------------------------------------------
# Copyright (c) 2018, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------
import numpy as np

from yt.frontends.stream.api import \
    load_uniform_grid


def test_point_periodicity():
    """Test for #1711."""
    domain_dims = (16, 16, 16)
    dens = np.ones(domain_dims)
    dens[3, 3, 3] = 100.0
    bbox = np.array([[-1.0, 1.0], [-1.0, 1.0], [-1.0, 1.0]])
    ds = load_uniform_grid({'density': dens}, domain_dims, length_unit='cm',
                           bbox=bbox, periodicity=(True, True, True))
    ad = ds.all_data()
    assert ds.r[ad.argmax('density')]['density'] == ad.max('density')
