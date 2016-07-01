from yt.testing import \
        fake_amr_ds, \
        assert_array_equal
import numpy as np

def get_ds():
    from yt.utilities.lib.geometry_utils import compute_morton
    def _morton_index(field, data):
        eps = np.finfo("f8").eps
        uq = data.ds.domain_left_edge.uq
        LE = data.ds.domain_left_edge - eps * uq
        RE = data.ds.domain_right_edge + eps * uq
        # .ravel() only copies if it needs to
        morton = compute_morton(data["index", "x"].ravel(),
                                data["index", "y"].ravel(),
                                data["index", "z"].ravel(), LE, RE)
        morton.shape = data["index", "x"].shape
        return morton.view("f8")
    ds = fake_amr_ds()
    ds.add_field(("index", "morton_index"), function=_morton_index,
                       units = "")
    return ds

def test_boolean_spheres_no_overlap():
    r"""Test to make sure that boolean objects (spheres, no overlap)
    behave the way we expect.

    Test non-overlapping spheres. This also checks that the original spheres
    don't change as part of constructing the booleans.
    """
    ds = get_ds()
    sp1 = ds.sphere([0.25, 0.25, 0.25], 0.15)
    sp2 = ds.sphere([0.75, 0.75, 0.75], 0.15)
    # Store the original indices
    i1 = sp1["index","morton_index"]
    i1.sort()
    i2 = sp2["index","morton_index"]
    i2.sort()
    ii = np.concatenate((i1, i2))
    ii.sort()
    # Make some booleans
    bo1 = sp1 & sp2
    bo2 = sp1 - sp2
    bo3 = sp1 | sp2 # also works with +
    # This makes sure the original containers didn't change.
    new_i1 = sp1["index","morton_index"]
    new_i1.sort()
    new_i2 = sp2["index","morton_index"]
    new_i2.sort()
    assert_array_equal(new_i1, i1)
    assert_array_equal(new_i2, i2)
    # Now make sure the indices also behave as we expect.
    empty = np.array([])
    assert_array_equal(bo1["index","morton_index"], empty)
    b2 = bo2["index","morton_index"]
    b2.sort()
    assert_array_equal(b2, i1)
    b3 = bo3["index","morton_index"]
    b3.sort()
    assert_array_equal(b3, ii)
 
def test_boolean_spheres_overlap():
    r"""Test to make sure that boolean objects (spheres, overlap)
    behave the way we expect.

    Test overlapping spheres.
    """
    ds = get_ds()
    sp1 = ds.sphere([0.45, 0.45, 0.45], 0.15)
    sp2 = ds.sphere([0.55, 0.55, 0.55], 0.15)
    # Get indices of both.
    i1 = sp1["index","morton_index"]
    i2 = sp2["index","morton_index"]
    # Make some booleans
    bo1 = sp1 & sp2
    bo2 = sp1 - sp2
    bo3 = sp1 | sp2
    # Now make sure the indices also behave as we expect.
    lens = np.intersect1d(i1, i2)
    apple = np.setdiff1d(i1, i2)
    both = np.union1d(i1, i2)
    b1 = bo1["index","morton_index"]
    b1.sort()
    b2 = bo2["index","morton_index"]
    b2.sort()
    b3 = bo3["index","morton_index"]
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
    ds = get_ds()
    re1 = ds.region([0.25]*3, [0.2]*3, [0.3]*3)
    re2 = ds.region([0.65]*3, [0.6]*3, [0.7]*3)
    # Store the original indices
    i1 = re1["index","morton_index"]
    i1.sort()
    i2 = re2["index","morton_index"]
    i2.sort()
    ii = np.concatenate((i1, i2))
    ii.sort()
    # Make some booleans
    bo1 = re1 & re2
    bo2 = re1 - re2
    bo3 = re1 | re2
    # This makes sure the original containers didn't change.
    new_i1 = re1["index","morton_index"]
    new_i1.sort()
    new_i2 = re2["index","morton_index"]
    new_i2.sort()
    assert_array_equal(new_i1, i1)
    assert_array_equal(new_i2, i2)
    # Now make sure the indices also behave as we expect.
    empty = np.array([])
    assert_array_equal(bo1["index","morton_index"], empty)
    b2 = bo2["index","morton_index"]
    b2.sort()
    assert_array_equal(b2, i1 )
    b3 = bo3["index","morton_index"]
    b3.sort()
    assert_array_equal(b3, ii)

def test_boolean_regions_overlap():
    r"""Test to make sure that boolean objects (regions, overlap)
    behave the way we expect.

    Test overlapping regions.
    """
    ds = get_ds()
    re1 = ds.region([0.55]*3, [0.5]*3, [0.6]*3)
    re2 = ds.region([0.6]*3, [0.55]*3, [0.65]*3)
    # Get indices of both.
    i1 = re1["index","morton_index"]
    i2 = re2["index","morton_index"]
    # Make some booleans
    bo1 = re1 & re2
    bo2 = re1 - re2
    bo3 = re1 | re2
    # Now make sure the indices also behave as we expect.
    cube = np.intersect1d(i1, i2)
    bite_cube = np.setdiff1d(i1, i2)
    both = np.union1d(i1, i2)
    b1 = bo1["index","morton_index"]
    b1.sort()
    b2 = bo2["index","morton_index"]
    b2.sort()
    b3 = bo3["index","morton_index"]
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
    ds = get_ds()
    cyl1 = ds.disk([0.25]*3, [1, 0, 0], 0.1, 0.1)
    cyl2 = ds.disk([0.75]*3, [1, 0, 0], 0.1, 0.1)
    # Store the original indices
    i1 = cyl1["index","morton_index"]
    i1.sort()
    i2 = cyl2["index","morton_index"]
    i2.sort()
    ii = np.concatenate((i1, i2))
    ii.sort()
    # Make some booleans
    bo1 = cyl1 & cyl2
    bo2 = cyl1 - cyl2
    bo3 = cyl1 | cyl2
    # This makes sure the original containers didn't change.
    new_i1 = cyl1["index","morton_index"]
    new_i1.sort()
    new_i2 = cyl2["index","morton_index"]
    new_i2.sort()
    assert_array_equal(new_i1, i1)
    assert_array_equal(new_i2, i2)
    # Now make sure the indices also behave as we expect.
    empty = np.array([])
    assert_array_equal(bo1["index","morton_index"], empty)
    b2 = bo2["index","morton_index"]
    b2.sort()
    assert_array_equal(b2, i1)
    b3 = bo3["index","morton_index"]
    b3.sort()
    assert_array_equal(b3, ii)

def test_boolean_cylinders_overlap():
    r"""Test to make sure that boolean objects (cylinders, overlap)
    behave the way we expect.

    Test overlapping cylinders.
    """
    ds = get_ds()
    cyl1 = ds.disk([0.45]*3, [1, 0, 0], 0.2, 0.2)
    cyl2 = ds.disk([0.55]*3, [1, 0, 0], 0.2, 0.2)
    # Get indices of both.
    i1 = cyl1["index","morton_index"]
    i2 = cyl2["index","morton_index"]
    # Make some booleans
    bo1 = cyl1 & cyl2
    bo2 = cyl1 - cyl2
    bo3 = cyl1 | cyl2
    # Now make sure the indices also behave as we expect.
    vlens = np.intersect1d(i1, i2)
    bite_disk = np.setdiff1d(i1, i2)
    both = np.union1d(i1, i2)
    b1 = bo1["index","morton_index"]
    b1.sort()
    b2 = bo2["index","morton_index"]
    b2.sort()
    b3 = bo3["index","morton_index"]
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
    ds = get_ds()
    ell1 = ds.ellipsoid([0.25]*3, 0.05, 0.05, 0.05, np.array([0.1]*3), 0.1)
    ell2 = ds.ellipsoid([0.75]*3, 0.05, 0.05, 0.05, np.array([0.1]*3), 0.1)
    # Store the original indices
    i1 = ell1["index","morton_index"]
    i1.sort()
    i2 = ell2["index","morton_index"]
    i2.sort()
    ii = np.concatenate((i1, i2))
    ii.sort()
    # Make some booleans
    bo1 = ell1 & ell2
    bo2 = ell1 - ell2
    bo3 = ell1 | ell2
    # This makes sure the original containers didn't change.
    new_i1 = ell1["index","morton_index"]
    new_i1.sort()
    new_i2 = ell2["index","morton_index"]
    new_i2.sort()
    assert_array_equal(new_i1, i1 )
    assert_array_equal(new_i2, i2)
    # Now make sure the indices also behave as we expect.
    empty = np.array([])
    assert_array_equal(bo1["index","morton_index"], empty)
    b2 = bo2["index","morton_index"]
    b2.sort()
    assert_array_equal(b2, i1)
    b3 = bo3["index","morton_index"]
    b3.sort()
    assert_array_equal(b3, ii)

def test_boolean_ellipsoids_overlap():
    r"""Test to make sure that boolean objects (ellipsoids, overlap)
    behave the way we expect.

    Test overlapping ellipsoids.
    """
    ds = get_ds()
    ell1 = ds.ellipsoid([0.45]*3, 0.05, 0.05, 0.05, np.array([0.1]*3), 0.1)
    ell2 = ds.ellipsoid([0.55]*3, 0.05, 0.05, 0.05, np.array([0.1]*3), 0.1)
    # Get indices of both.
    i1 = ell1["index","morton_index"]
    i2 = ell2["index","morton_index"]
    # Make some booleans
    bo1 = ell1 & ell2
    bo2 = ell1 - ell2
    bo3 = ell1 | ell2
    # Now make sure the indices also behave as we expect.
    overlap = np.intersect1d(i1, i2)
    diff = np.setdiff1d(i1, i2)
    both = np.union1d(i1, i2)
    b1 = bo1["index","morton_index"]
    b1.sort()
    b2 = bo2["index","morton_index"]
    b2.sort()
    b3 = bo3["index","morton_index"]
    b3.sort()
    assert_array_equal(b1, overlap)
    assert_array_equal(b2, diff)
    assert_array_equal(b3, both)

def test_boolean_mix_periodicity():
    r"""Test that a hybrid boolean region behaves as we expect.

    This also tests nested logic and that periodicity works.
    """
    ds = get_ds()
    re = ds.region([0.5]*3, [0.0]*3, [1]*3) # whole thing
    sp = ds.sphere([0.95]*3, 0.3) # wraps around
    cyl = ds.disk([0.05]*3, [1,1,1], 0.1, 0.4) # wraps around
    # Get original indices
    rei = re["index","morton_index"]
    spi = sp["index","morton_index"]
    cyli = cyl["index","morton_index"]
    # Make some booleans
    # whole box minux spherical bites at corners
    bo1 = re - sp
    # sphere plus cylinder
    bo2 = sp | cyl
    # a jumble, the region minus the sp+cyl
    bo3 = re - (sp | cyl)
    # Now make sure the indices also behave as we expect.
    expect = np.setdiff1d(rei, spi)
    ii = bo1["index","morton_index"]
    ii.sort()
    assert_array_equal(expect, ii)
    #
    expect = np.union1d(spi, cyli)
    ii = bo2["index","morton_index"]
    ii.sort()
    assert_array_equal(expect, ii)
    #
    expect = np.union1d(spi, cyli)
    expect = np.setdiff1d(rei, expect)
    ii = bo3["index","morton_index"]
    ii.sort()
    assert_array_equal(expect, ii)

