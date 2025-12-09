# Some tests for finding bounding boxes

import numpy as np
from numpy.testing import assert_equal

from yt.testing import assert_allclose_units, fake_amr_ds, fake_octree_ds


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


def test_covering_grid_bbox():
    ds = fake_octree_ds(num_zones=32)

    cg = ds.covering_grid(level=0, left_edge=[0.3, 0.3, 0.3], dims=[8, 8, 8])

    # Make a covering grid with a data source
    sp = ds.sphere(
        [0.5, 0.5, 0.5],
        0.15,
    )
    cgds = ds.covering_grid(
        level=0,
        left_edge=[0.3, 0.3, 0.3],
        dims=[8, 8, 8],
        data_source=sp,
    )

    # The bounding boxes of the two covering grids should be the same
    cg_bbox = cg.get_bbox()
    cgds_bbox = cgds.get_bbox()
    assert_equal(cg_bbox, cgds_bbox)

    # The bounding box of the cg's data source should be the left edge of sp and right edge of cgds
    cgds_ds_bbox = cgds._data_source.get_bbox()
    le_sp, re_sp = sp.get_bbox()

    assert_equal(le_sp, cgds_ds_bbox[0])
    assert_equal(cgds_bbox[1], cgds_ds_bbox[1])


def test_intersection_bbox():
    ds = fake_octree_ds(num_zones=32)

    # intersect a region and a sphere smaller than the region, get back
    # the bbox of the sphere
    sp1 = ds.sphere(ds.domain_center, (0.1, "unitary"))
    reg = ds.region(ds.domain_center, ds.domain_left_edge, ds.domain_right_edge)

    le, re = ds.intersection((sp1, reg)).get_bbox()
    le_sp, re_sp = sp1.get_bbox()
    assert_equal(le_sp, le)
    assert_equal(re_sp, re)

    # check for no error with a single data source
    le, re = ds.intersection((reg,)).get_bbox()
    assert_equal(le, ds.domain_left_edge)
    assert_equal(re, ds.domain_right_edge)

    # check some overlapping regions shifted along one dimension
    c = ds.domain_center - ds.arr((0.3, 0, 0), "code_length")
    le = c - ds.quan(0.2, "code_length")
    re = c + ds.quan(0.2, "code_length")
    r1 = ds.region(c, le, re)
    offset = ds.arr((0.1, 0, 0), "code_length")
    r2 = ds.region(c + offset, le + offset, re + offset)
    r3 = ds.region(c + 2 * offset, le + 2 * offset, re + 2 * offset)

    r_int = ds.intersection((r1, r2, r3))
    le, re = r_int.get_bbox()
    assert_equal(le, ds.arr((0.2, 0.3, 0.3), "code_length"))
    assert_equal(re, ds.arr((0.4, 0.7, 0.7), "code_length"))
