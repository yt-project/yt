from yt.testing import *
from yt.data_objects.api import add_field

def setup():
    from yt.config import ytcfg
    ytcfg["yt","__withintesting"] = "True"
    def _ID(field, data):
        width = data.pf.domain_right_edge - data.pf.domain_left_edge
        delta = width / data.pf.h.get_smallest_dx()
        x = data['x'] - data.pf.h.get_smallest_dx() / 2.
        y = data['y'] - data.pf.h.get_smallest_dx() / 2.
        z = data['z'] - data.pf.h.get_smallest_dx() / 2.
        xi = x / data.pf.h.get_smallest_dx()
        yi = y / data.pf.h.get_smallest_dx()
        zi = z / data.pf.h.get_smallest_dx()
        index = xi + delta[0] * (yi + delta[1] * zi)
        index = index.astype('int64')
        return index

    add_field("ID", function=_ID)

def test_boolean_spheres_no_overlap():
    r"""Test to make sure that boolean objects (spheres, no overlap)
    behave the way we expect.

    Test non-overlapping spheres. This also checks that the original spheres
    don't change as part of constructing the booleans.
    """
    pf = fake_random_pf(64)
    pf.h
    sp1 = pf.h.sphere([0.25, 0.25, 0.25], 0.15)
    sp2 = pf.h.sphere([0.75, 0.75, 0.75], 0.15)
    # Store the original indices
    i1 = sp1['ID']
    i1.sort()
    i2 = sp2['ID']
    i2.sort()
    ii = np.concatenate((i1, i2))
    ii.sort()
    # Make some booleans
    bo1 = pf.h.boolean([sp1, "AND", sp2]) # empty
    bo2 = pf.h.boolean([sp1, "NOT", sp2]) # only sp1
    bo3 = pf.h.boolean([sp1, "OR", sp2]) # combination
    # This makes sure the original containers didn't change.
    new_i1 = sp1['ID']
    new_i1.sort()
    new_i2 = sp2['ID']
    new_i2.sort()
    assert_equal(new_i1, i1)
    assert_equal(new_i2, i2)
    # Now make sure the indices also behave as we expect.
    empty = np.array([])
    assert_array_equal(bo1['ID'], empty)
    b2 = bo2['ID']
    b2.sort()
    assert_array_equal(b2, i1)
    b3 = bo3['ID']
    b3.sort()
    assert_array_equal(b3, ii)
 
def test_boolean_spheres_overlap():
    r"""Test to make sure that boolean objects (spheres, overlap)
    behave the way we expect.

    Test overlapping spheres.
    """
    pf = fake_random_pf(64)
    pf.h
    sp1 = pf.h.sphere([0.45, 0.45, 0.45], 0.15)
    sp2 = pf.h.sphere([0.55, 0.55, 0.55], 0.15)
    # Get indices of both.
    i1 = sp1['ID']
    i2 = sp2['ID']
    # Make some booleans
    bo1 = pf.h.boolean([sp1, "AND", sp2]) # overlap (a lens)
    bo2 = pf.h.boolean([sp1, "NOT", sp2]) # sp1 - sp2 (sphere with bite)
    bo3 = pf.h.boolean([sp1, "OR", sp2]) # combination (H2)
    # Now make sure the indices also behave as we expect.
    lens = np.intersect1d(i1, i2)
    apple = np.setdiff1d(i1, i2)
    both = np.union1d(i1, i2)
    b1 = bo1['ID']
    b1.sort()
    b2 = bo2['ID']
    b2.sort()
    b3 = bo3['ID']
    b3.sort()
    assert_array_equal(b1, lens)
    assert_array_equal(b2, apple)
    assert_array_equal(b3, both)

def test_boolean_regions_no_overlap():
    r"""Test to make sure that boolean objects (regions, no overlap)
    behave the way we expect.

    Test non-overlapping regions. This also checks that the original regions
    don't change as part of constructing the booleans.
    """
    pf = fake_random_pf(64)
    pf.h
    re1 = pf.h.region([0.25]*3, [0.2]*3, [0.3]*3)
    re2 = pf.h.region([0.65]*3, [0.6]*3, [0.7]*3)
    # Store the original indices
    i1 = re1['ID']
    i1.sort()
    i2 = re2['ID']
    i2.sort()
    ii = np.concatenate((i1, i2))
    ii.sort()
    # Make some booleans
    bo1 = pf.h.boolean([re1, "AND", re2]) # empty
    bo2 = pf.h.boolean([re1, "NOT", re2]) # only re1
    bo3 = pf.h.boolean([re1, "OR", re2]) # combination
    # This makes sure the original containers didn't change.
    new_i1 = re1['ID']
    new_i1.sort()
    new_i2 = re2['ID']
    new_i2.sort()
    assert_equal(new_i1, i1)
    assert_equal(new_i2, i2)
    # Now make sure the indices also behave as we expect.
    empty = np.array([])
    assert_array_equal(bo1['ID'], empty)
    b2 = bo2['ID']
    b2.sort()
    assert_array_equal(b2, i1)
    b3 = bo3['ID']
    b3.sort()
    assert_array_equal(b3, ii)

def test_boolean_regions_overlap():
    r"""Test to make sure that boolean objects (regions, overlap)
    behave the way we expect.

    Test overlapping regions.
    """
    pf = fake_random_pf(64)
    pf.h
    re1 = pf.h.region([0.55]*3, [0.5]*3, [0.6]*3)
    re2 = pf.h.region([0.6]*3, [0.55]*3, [0.65]*3)
    # Get indices of both.
    i1 = re1['ID']
    i2 = re2['ID']
    # Make some booleans
    bo1 = pf.h.boolean([re1, "AND", re2]) # overlap (small cube)
    bo2 = pf.h.boolean([re1, "NOT", re2]) # sp1 - sp2 (large cube with bite)
    bo3 = pf.h.boolean([re1, "OR", re2]) # combination (merged large cubes)
    # Now make sure the indices also behave as we expect.
    cube = np.intersect1d(i1, i2)
    bite_cube = np.setdiff1d(i1, i2)
    both = np.union1d(i1, i2)
    b1 = bo1['ID']
    b1.sort()
    b2 = bo2['ID']
    b2.sort()
    b3 = bo3['ID']
    b3.sort()
    assert_array_equal(b1, cube)
    assert_array_equal(b2, bite_cube)
    assert_array_equal(b3, both)

def test_boolean_cylinders_no_overlap():
    r"""Test to make sure that boolean objects (cylinders, no overlap)
    behave the way we expect.

    Test non-overlapping cylinders. This also checks that the original cylinders
    don't change as part of constructing the booleans.
    """
    pf = fake_random_pf(64)
    pf.h
    cyl1 = pf.h.disk([0.25]*3, [1, 0, 0], 0.1, 0.1)
    cyl2 = pf.h.disk([0.75]*3, [1, 0, 0], 0.1, 0.1)
    # Store the original indices
    i1 = cyl1['ID']
    i1.sort()
    i2 = cyl2['ID']
    i2.sort()
    ii = np.concatenate((i1, i2))
    ii.sort()
    # Make some booleans
    bo1 = pf.h.boolean([cyl1, "AND", cyl2]) # empty
    bo2 = pf.h.boolean([cyl1, "NOT", cyl2]) # only cyl1
    bo3 = pf.h.boolean([cyl1, "OR", cyl2]) # combination
    # This makes sure the original containers didn't change.
    new_i1 = cyl1['ID']
    new_i1.sort()
    new_i2 = cyl2['ID']
    new_i2.sort()
    assert_equal(new_i1, i1)
    assert_equal(new_i2, i2)
    # Now make sure the indices also behave as we expect.
    empty = np.array([])
    assert_array_equal(bo1['ID'], empty)
    b2 = bo2['ID']
    b2.sort()
    assert_array_equal(b2, i1)
    b3 = bo3['ID']
    b3.sort()
    assert_array_equal(b3, ii)

def test_boolean_cylinders_overlap():
    r"""Test to make sure that boolean objects (cylinders, overlap)
    behave the way we expect.

    Test overlapping cylinders.
    """
    pf = fake_random_pf(64)
    pf.h
    cyl1 = pf.h.disk([0.45]*3, [1, 0, 0], 0.2, 0.2)
    cyl2 = pf.h.disk([0.55]*3, [1, 0, 0], 0.2, 0.2)
    # Get indices of both.
    i1 = cyl1['ID']
    i2 = cyl2['ID']
    # Make some booleans
    bo1 = pf.h.boolean([cyl1, "AND", cyl2]) # overlap (vertically extened lens)
    bo2 = pf.h.boolean([cyl1, "NOT", cyl2]) # sp1 - sp2 (disk minus a bite)
    bo3 = pf.h.boolean([cyl1, "OR", cyl2]) # combination (merged disks)
    # Now make sure the indices also behave as we expect.
    vlens = np.intersect1d(i1, i2)
    bite_disk = np.setdiff1d(i1, i2)
    both = np.union1d(i1, i2)
    b1 = bo1['ID']
    b1.sort()
    b2 = bo2['ID']
    b2.sort()
    b3 = bo3['ID']
    b3.sort()
    assert_array_equal(b1, vlens)
    assert_array_equal(b2, bite_disk)
    assert_array_equal(b3, both)

def test_boolean_ellipsoids_no_overlap():
    r"""Test to make sure that boolean objects (ellipsoids, no overlap)
    behave the way we expect.

    Test non-overlapping ellipsoids. This also checks that the original
    ellipsoids don't change as part of constructing the booleans.
    """
    pf = fake_random_pf(64)
    pf.h
    ell1 = pf.h.ellipsoid([0.25]*3, 0.05, 0.05, 0.05, np.array([0.1]*3),
        np.array([0.1]*3))
    ell2 = pf.h.ellipsoid([0.75]*3, 0.05, 0.05, 0.05, np.array([0.1]*3),
        np.array([0.1]*3))
    # Store the original indices
    i1 = ell1['ID']
    i1.sort()
    i2 = ell2['ID']
    i2.sort()
    ii = np.concatenate((i1, i2))
    ii.sort()
    # Make some booleans
    bo1 = pf.h.boolean([ell1, "AND", ell2]) # empty
    bo2 = pf.h.boolean([ell1, "NOT", ell2]) # only cyl1
    bo3 = pf.h.boolean([ell1, "OR", ell2]) # combination
    # This makes sure the original containers didn't change.
    new_i1 = ell1['ID']
    new_i1.sort()
    new_i2 = ell2['ID']
    new_i2.sort()
    assert_equal(new_i1, i1)
    assert_equal(new_i2, i2)
    # Now make sure the indices also behave as we expect.
    empty = np.array([])
    assert_array_equal(bo1['ID'], empty)
    b2 = bo2['ID']
    b2.sort()
    assert_array_equal(b2, i1)
    b3 = bo3['ID']
    b3.sort()
    assert_array_equal(b3, ii)

def test_boolean_ellipsoids_overlap():
    r"""Test to make sure that boolean objects (ellipsoids, overlap)
    behave the way we expect.

    Test overlapping ellipsoids.
    """
    pf = fake_random_pf(64)
    pf.h
    ell1 = pf.h.ellipsoid([0.45]*3, 0.05, 0.05, 0.05, np.array([0.1]*3),
        np.array([0.1]*3))
    ell2 = pf.h.ellipsoid([0.55]*3, 0.05, 0.05, 0.05, np.array([0.1]*3),
        np.array([0.1]*3))
    # Get indices of both.
    i1 = ell1['ID']
    i2 = ell2['ID']
    # Make some booleans
    bo1 = pf.h.boolean([ell1, "AND", ell2]) # overlap
    bo2 = pf.h.boolean([ell1, "NOT", ell2]) # ell1 - ell2
    bo3 = pf.h.boolean([ell1, "OR", ell2]) # combination
    # Now make sure the indices also behave as we expect.
    overlap = np.intersect1d(i1, i2)
    diff = np.setdiff1d(i1, i2)
    both = np.union1d(i1, i2)
    b1 = bo1['ID']
    b1.sort()
    b2 = bo2['ID']
    b2.sort()
    b3 = bo3['ID']
    b3.sort()
    assert_array_equal(b1, overlap)
    assert_array_equal(b2, diff)
    assert_array_equal(b3, both)

def test_boolean_mix_periodicity():
    r"""Test that a hybrid boolean region behaves as we expect.

    This also tests nested logic and that periodicity works.
    """
    pf = fake_random_pf(64)
    pf.h
    re = pf.h.region([0.5]*3, [0.0]*3, [1]*3) # whole thing
    sp = pf.h.sphere([0.95]*3, 0.3) # wraps around
    cyl = pf.h.disk([0.05]*3, [1,1,1], 0.1, 0.4) # wraps around
    ell = pf.h.ellipsoid([0.35]*3, 0.05, 0.05, 0.05, np.array([0.1]*3),
        np.array([0.1]*3)) # no wrap
    # Get original indices
    rei = re['ID']
    spi = sp['ID']
    cyli = cyl['ID']
    elli = ell['ID']
    # Make some booleans
    # whole box minux spherical bites at corners
    bo1 = pf.h.boolean([re, "NOT", sp])
    # sphere plus cylinder
    bo2 = pf.h.boolean([sp, "OR", cyl])
    # a big jumble, the region minus the ell+cyl (scepter shaped?), plus the
    # sphere which should add back some of what the ell+cyl took out.
    bo3 = pf.h.boolean([re, "NOT", "(", ell, "OR", cyl, ")", "OR", sp])
    # Now make sure the indices also behave as we expect.
    expect = np.setdiff1d(rei, spi)
    ii = bo1['ID']
    ii.sort()
    assert_array_equal(expect, ii)
    #
    expect = np.union1d(spi, cyli)
    ii = bo2['ID']
    ii.sort()
    assert_array_equal(expect, ii)
    #
    expect = np.union1d(elli, cyli)
    expect = np.setdiff1d(rei, expect)
    expect = np.union1d(expect, spi)
    ii = bo3['ID']
    ii.sort()
    assert_array_equal(expect, ii)

