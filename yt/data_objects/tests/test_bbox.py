# Some tests for finding bounding boxes

import numpy as np

from yt.testing import \
    fake_amr_ds, \
    assert_equal, \
    assert_allclose_units


def test_object_bbox():
    ds = fake_amr_ds()
    reg = ds.box(ds.domain_left_edge+0.5*ds.domain_width,
                 ds.domain_right_edge-0.5*ds.domain_width)
    le, re = reg.get_bbox()
    assert_equal(le, ds.domain_left_edge+0.5*ds.domain_width)
    assert_equal(re, ds.domain_right_edge-0.5*ds.domain_width)
    sp = ds.sphere("c", (0.1, "unitary"))
    le, re = sp.get_bbox()
    assert_equal(le, -sp.radius+sp.center)
    assert_equal(re, sp.radius+sp.center)
    dk = ds.disk("c", [1,1,0], (0.25, "unitary"), (0.25, "unitary"))
    le, re = dk.get_bbox()
    le0 = ds.arr([0.5-0.25*np.sqrt(2.0), 0.5-0.25*np.sqrt(2.0), 0.25], "code_length")
    re0 = ds.arr([0.5+0.25*np.sqrt(2.0), 0.5+0.25*np.sqrt(2.0), 0.75], "code_length")
    assert_allclose_units(le, le0)
    assert_allclose_units(re, re0)
    ep = ds.ellipsoid("c", 0.3, 0.2, 0.1, np.array([0.1, 0.1, 0.1]), 0.2)
    le, re = ep.get_bbox()
    assert_equal(le, -ds.quan(0.3, "code_length")+sp.center)
    assert_equal(re, ds.quan(0.3, "code_length")+sp.center)

