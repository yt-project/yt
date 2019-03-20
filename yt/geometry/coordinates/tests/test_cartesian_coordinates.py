# Some tests for the Cartesian coordinates handler

import numpy as np

from yt.testing import \
    fake_amr_ds, \
    assert_equal, \
    assert_allclose_units

# Our canonical tests are that we can access all of our fields and we can
# compute our volume correctly.

def test_cartesian_coordinates():
    # We're going to load up a simple AMR grid and check its volume
    # calculations and path length calculations.
    ds = fake_amr_ds()
    axes = list(set(ds.coordinates.axis_name.values()))
    axes.sort()
    for i, axis in enumerate(axes):
        dd = ds.all_data()
        fi = ("index", axis)
        fd = ("index", "d%s" % axis)
        fp = ("index", "path_element_%s" % axis)
        ma = np.argmax(dd[fi])
        assert_equal(dd[fi][ma] + dd[fd][ma] / 2.0, ds.domain_right_edge[i])
        mi = np.argmin(dd[fi])
        assert_equal(dd[fi][mi] - dd[fd][mi] / 2.0, ds.domain_left_edge[i])
        assert_equal(dd[fd].min(), ds.index.get_smallest_dx())
        assert_equal(dd[fd].max(), (ds.domain_width/ds.domain_dimensions)[i])
        assert_equal(dd[fd], dd[fp])
    assert_equal(dd["cell_volume"].sum(dtype="float64"), ds.domain_width.prod())


def test_object_bbox():
    from yt.geometry.coordinates.cartesian_coordinates \
        import _find_object_bbox
    ds = fake_amr_ds()
    reg = ds.box(ds.domain_left_edge+0.5*ds.domain_width,
                 ds.domain_right_edge-0.5*ds.domain_width)
    le, re = _find_object_bbox(reg)
    assert_equal(le, ds.domain_left_edge+0.5*ds.domain_width)
    assert_equal(re, ds.domain_right_edge-0.5*ds.domain_width)
    sp = ds.sphere("c", (0.1, "unitary"))
    le, re = _find_object_bbox(sp)
    assert_equal(le, -sp.radius+sp.center)
    assert_equal(re, sp.radius+sp.center)
    ray = ds.ray(ds.domain_left_edge+0.5*ds.domain_width,
                 ds.domain_right_edge-0.5*ds.domain_width)
    le, re = _find_object_bbox(ray)
    assert_equal(le, ds.domain_left_edge+0.5*ds.domain_width)
    assert_equal(re, ds.domain_right_edge-0.5*ds.domain_width)
    oray = ds.ortho_ray(1, (0.2*ds.domain_width[0]+ds.domain_left_edge[0], 
                            0.7*ds.domain_width[2]+ds.domain_left_edge[0]))
    le, re = _find_object_bbox(oray)
    le0 = ds.domain_left_edge.copy()
    re0 = ds.domain_right_edge.copy()
    le0[0] = 0.7*ds.domain_width[0]+ds.domain_left_edge[0]
    le0[2] = 0.2*ds.domain_width[2]+ds.domain_left_edge[2]
    re0[0] = 0.7*ds.domain_width[0]+ds.domain_left_edge[0]
    re0[2] = 0.2*ds.domain_width[2]+ds.domain_left_edge[2]
    assert_equal(le, le0)
    assert_equal(re, re0)
    dk = ds.disk("c", [1,1,0], (0.25, "unitary"), (0.25, "unitary"))
    le, re = _find_object_bbox(dk)
    le0 = ds.arr([0.5-0.25*np.sqrt(2.0), 0.5-0.25*np.sqrt(2.0), 0.25], "code_length")
    re0 = ds.arr([0.5+0.25*np.sqrt(2.0), 0.5+0.25*np.sqrt(2.0), 0.75], "code_length")
    assert_allclose_units(le, le0)
    assert_allclose_units(re, re0)
    ep = ds.ellipsoid("c", 0.3, 0.2, 0.1, np.array([0.1, 0.1, 0.1]), 0.2)
    le, re = _find_object_bbox(ep)
    assert_equal(le, -ds.quan(0.3, "code_length")+sp.center)
    assert_equal(re, ds.quan(0.3, "code_length")+sp.center)

