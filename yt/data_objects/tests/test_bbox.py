# Some tests for finding bounding boxes

import numpy as np

from yt.testing import assert_allclose_units, assert_equal, fake_amr_ds


def test_object_bbox():
    ds = fake_amr_ds()
    reg = ds.box(
        ds.domain_left_edge + 0.5 * ds.domain_width,
        ds.domain_right_edge - 0.5 * ds.domain_width,
    )
    le, re = reg.get_bbox()
    assert_equal(le, ds.domain_left_edge + 0.5 * ds.domain_width)
    assert_equal(re, ds.domain_right_edge - 0.5 * ds.domain_width)
    sp = ds.sphere("c", (0.1, "unitary"))
    le, re = sp.get_bbox()
    assert_equal(le, -sp.radius + sp.center)
    assert_equal(re, sp.radius + sp.center)
    dk = ds.disk("c", [1, 1, 0], (0.25, "unitary"), (0.25, "unitary"))
    le, re = dk.get_bbox()
    le0 = ds.arr(
        [0.5 - 0.25 * np.sqrt(2.0), 0.5 - 0.25 * np.sqrt(2.0), 0.25], "code_length"
    )
    re0 = ds.arr(
        [0.5 + 0.25 * np.sqrt(2.0), 0.5 + 0.25 * np.sqrt(2.0), 0.75], "code_length"
    )
    assert_allclose_units(le, le0)
    assert_allclose_units(re, re0)
    ep = ds.ellipsoid("c", 0.3, 0.2, 0.1, np.array([0.1, 0.1, 0.1]), 0.2)
    le, re = ep.get_bbox()
    assert_equal(le, -ds.quan(0.3, "code_length") + sp.center)
    assert_equal(re, ds.quan(0.3, "code_length") + sp.center)
    spb = ds.sphere(
        ds.domain_center - ds.quan(0.1, "code_length"), (0.1, "code_length")
    )
    regb = ds.box(ds.domain_center, ds.domain_center + ds.quan(0.2, "code_length"))
    br1 = spb & regb
    br2 = spb | regb
    br3 = spb ^ regb
    br4 = ~regb
    le1, re1 = br1.get_bbox()
    le2, re2 = br2.get_bbox()
    le3, re3 = br3.get_bbox()
    le4, re4 = br4.get_bbox()
    le0 = ds.arr([0.3, 0.3, 0.3], "code_length")
    re0 = ds.arr([0.7, 0.7, 0.7], "code_length")
    assert_allclose_units(le1, le0)
    assert_allclose_units(re1, re0)
    assert_allclose_units(le2, le0)
    assert_allclose_units(re2, re0)
    assert_allclose_units(le3, le0)
    assert_allclose_units(re3, re0)
    assert_equal(le4, regb.left_edge)
    assert_equal(re4, regb.right_edge)
