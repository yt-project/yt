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
from yt.funcs import ensure_tuple
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
    for resolution in [64, (64, 64), (32, 64), (64, 32)]:
        frb = slc.to_frb((0.5,"unitary"), resolution)
        frb_ds = frb.export_dataset(fields=["density"], nprocs=8)
        dd_frb = frb_ds.all_data()

        assert_equal(frb_ds.domain_left_edge.v, np.array([0.25,0.25,0.0]))
        assert_equal(frb_ds.domain_right_edge.v, np.array([0.75,0.75,1.0]))
        assert_equal(frb_ds.domain_width.v, np.array([0.5,0.5,1.0]))
        if len(ensure_tuple(resolution)) == 1:
            resolution = (resolution, resolution)
        frb_ds_dd = np.array([resolution[0], resolution[1], 1], dtype="int64")
        assert_equal(frb_ds.domain_dimensions, frb_ds_dd)
        dens_image = frb["density"]
        assert_allclose_units(dens_image.sum(),
                              dd_frb.quantities.total_quantity("density"))
        assert_equal(dens_image.shape, resolution)
        assert_equal(frb_ds.index.num_grids, 8)
