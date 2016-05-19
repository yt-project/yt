"""
Tests for exporting an FRB as a dataset



"""
from __future__ import absolute_import

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------
import numpy as np
from yt.testing import \
    fake_random_ds, assert_equal, \
    assert_allclose_units

def setup():
    """Test specific setup."""
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"


def test_export_frb():
    test_ds = fake_random_ds(128)
    slc = test_ds.slice(0,0.5)
    frb = slc.to_frb((0.5,"unitary"), 64)
    frb_ds = frb.export_dataset(fields=["density"], nprocs=8)
    dd_frb = frb_ds.all_data()

    yield assert_equal, frb_ds.domain_left_edge.v, np.array([0.25,0.25,0.0])
    yield assert_equal, frb_ds.domain_right_edge.v, np.array([0.75,0.75,1.0])
    yield assert_equal, frb_ds.domain_width.v, np.array([0.5,0.5,1.0])
    yield assert_equal, frb_ds.domain_dimensions, np.array([64,64,1], dtype="int64")
    yield assert_allclose_units, frb["density"].sum(), \
        dd_frb.quantities.total_quantity("density")
    yield assert_equal, frb_ds.index.num_grids, 8
