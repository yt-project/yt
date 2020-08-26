import numpy as np

from yt.testing import assert_array_equal, fake_amr_ds

# We use morton indices in this test because they are single floating point
# values that uniquely identify each cell.  That's a convenient way to compare
# inclusion in set operations, since there are no duplicates.


def test_boolean_spheres_no_overlap():
    r"""Test to make sure that boolean objects (spheres, no overlap)
    behave the way we expect.

    Test non-overlapping spheres. This also checks that the original spheres
    don't change as part of constructing the booleans.
    """
    ds = fake_amr_ds()
    sp1 = ds.sphere([0.25, 0.25, 0.25], 0.15)
    sp2 = ds.sphere([0.75, 0.75, 0.75], 0.15)
    # Store the original indices
    i1 = sp1["index", "morton_index"]
    i1.sort()
    i2 = sp2["index", "morton_index"]
    i2.sort()
    ii = np.concatenate((i1, i2))
    ii.sort()
    # Make some booleans
    bo1 = sp1 & sp2
    bo2 = sp1 - sp2
    bo3 = sp1 | sp2  # also works with +
    bo4 = ds.union([sp1, sp2])
    bo5 = ds.intersection([sp1, sp2])
    # This makes sure the original containers didn't change.
    new_i1 = sp1["index", "morton_index"]
    new_i1.sort()
    new_i2 = sp2["index", "morton_index"]
    new_i2.sort()
    assert_array_equal(new_i1, i1)
    assert_array_equal(new_i2, i2)
    # Now make sure the indices also behave as we expect.
    empty = np.array([])
    assert_array_equal(bo1["index", "morton_index"], empty)
    assert_array_equal(bo5["index", "morton_index"], empty)
    b2 = bo2["index", "morton_index"]
    b2.sort()
    assert_array_equal(b2, i1)
    b3 = bo3["index", "morton_index"]
    b3.sort()
    assert_array_equal(b3, ii)
    b4 = bo4["index", "morton_index"]
    b4.sort()
    assert_array_equal(b4, ii)
    bo6 = sp1 ^ sp2
    b6 = bo6["index", "morton_index"]
    b6.sort()
    assert_array_equal(b6, np.setxor1d(i1, i2))


def test_boolean_spheres_overlap():
    r"""Test to make sure that boolean objects (spheres, overlap)
    behave the way we expect.

    Test overlapping spheres.
    """
    ds = fake_amr_ds()
    sp1 = ds.sphere([0.45, 0.45, 0.45], 0.15)
    sp2 = ds.sphere([0.55, 0.55, 0.55], 0.15)
    # Get indices of both.
    i1 = sp1["index", "morton_index"]
    i2 = sp2["index", "morton_index"]
    # Make some booleans
    bo1 = sp1 & sp2
    bo2 = sp1 - sp2
    bo3 = sp1 | sp2
    bo4 = ds.union([sp1, sp2])
    bo5 = ds.intersection([sp1, sp2])
    # Now make sure the indices also behave as we expect.
    lens = np.intersect1d(i1, i2)
    apple = np.setdiff1d(i1, i2)
    both = np.union1d(i1, i2)
    b1 = bo1["index", "morton_index"]
    b1.sort()
    b2 = bo2["index", "morton_index"]
    b2.sort()
    b3 = bo3["index", "morton_index"]
    b3.sort()
    assert_array_equal(b1, lens)
    assert_array_equal(b2, apple)
    assert_array_equal(b3, both)
    b4 = bo4["index", "morton_index"]
    b4.sort()
    b5 = bo5["index", "morton_index"]
    b5.sort()
    assert_array_equal(b3, b4)
    assert_array_equal(b1, b5)
    bo6 = sp1 ^ sp2
    b6 = bo6["index", "morton_index"]
    b6.sort()
    assert_array_equal(b6, np.setxor1d(i1, i2))


def test_boolean_regions_no_overlap():
    r"""Test to make sure that boolean objects (regions, no overlap)
    behave the way we expect.

    Test non-overlapping regions. This also checks that the original regions
    don't change as part of constructing the booleans.
    """
    ds = fake_amr_ds()
    re1 = ds.region([0.25] * 3, [0.2] * 3, [0.3] * 3)
    re2 = ds.region([0.65] * 3, [0.6] * 3, [0.7] * 3)
    # Store the original indices
    i1 = re1["index", "morton_index"]
    i1.sort()
    i2 = re2["index", "morton_index"]
    i2.sort()
    ii = np.concatenate((i1, i2))
    ii.sort()
    # Make some booleans
    bo1 = re1 & re2
    bo2 = re1 - re2
    bo3 = re1 | re2
    bo4 = ds.union([re1, re2])
    bo5 = ds.intersection([re1, re2])
    # This makes sure the original containers didn't change.
    new_i1 = re1["index", "morton_index"]
    new_i1.sort()
    new_i2 = re2["index", "morton_index"]
    new_i2.sort()
    assert_array_equal(new_i1, i1)
    assert_array_equal(new_i2, i2)
    # Now make sure the indices also behave as we expect.
    empty = np.array([])
    assert_array_equal(bo1["index", "morton_index"], empty)
    assert_array_equal(bo5["index", "morton_index"], empty)
    b2 = bo2["index", "morton_index"]
    b2.sort()
    assert_array_equal(b2, i1)
    b3 = bo3["index", "morton_index"]
    b3.sort()
    assert_array_equal(b3, ii)
    b4 = bo4["index", "morton_index"]
    b4.sort()
    b5 = bo5["index", "morton_index"]
    b5.sort()
    assert_array_equal(b3, b4)
    bo6 = re1 ^ re2
    b6 = bo6["index", "morton_index"]
    b6.sort()
    assert_array_equal(b6, np.setxor1d(i1, i2))


def test_boolean_regions_overlap():
    r"""Test to make sure that boolean objects (regions, overlap)
    behave the way we expect.

    Test overlapping regions.
    """
    ds = fake_amr_ds()
    re1 = ds.region([0.55] * 3, [0.5] * 3, [0.6] * 3)
    re2 = ds.region([0.6] * 3, [0.55] * 3, [0.65] * 3)
    # Get indices of both.
    i1 = re1["index", "morton_index"]
    i2 = re2["index", "morton_index"]
    # Make some booleans
    bo1 = re1 & re2
    bo2 = re1 - re2
    bo3 = re1 | re2
    bo4 = ds.union([re1, re2])
    bo5 = ds.intersection([re1, re2])
    # Now make sure the indices also behave as we expect.
    cube = np.intersect1d(i1, i2)
    bite_cube = np.setdiff1d(i1, i2)
    both = np.union1d(i1, i2)
    b1 = bo1["index", "morton_index"]
    b1.sort()
    b2 = bo2["index", "morton_index"]
    b2.sort()
    b3 = bo3["index", "morton_index"]
    b3.sort()
    assert_array_equal(b1, cube)
    assert_array_equal(b2, bite_cube)
    assert_array_equal(b3, both)
    b4 = bo4["index", "morton_index"]
    b4.sort()
    b5 = bo5["index", "morton_index"]
    b5.sort()
    assert_array_equal(b3, b4)
    assert_array_equal(b1, b5)
    bo6 = re1 ^ re2
    b6 = bo6["index", "morton_index"]
    b6.sort()
    assert_array_equal(b6, np.setxor1d(i1, i2))


def test_boolean_cylinders_no_overlap():
    r"""Test to make sure that boolean objects (cylinders, no overlap)
    behave the way we expect.

    Test non-overlapping cylinders. This also checks that the original cylinders
    don't change as part of constructing the booleans.
    """
    ds = fake_amr_ds()
    cyl1 = ds.disk([0.25] * 3, [1, 0, 0], 0.1, 0.1)
    cyl2 = ds.disk([0.75] * 3, [1, 0, 0], 0.1, 0.1)
    # Store the original indices
    i1 = cyl1["index", "morton_index"]
    i1.sort()
    i2 = cyl2["index", "morton_index"]
    i2.sort()
    ii = np.concatenate((i1, i2))
    ii.sort()
    # Make some booleans
    bo1 = cyl1 & cyl2
    bo2 = cyl1 - cyl2
    bo3 = cyl1 | cyl2
    bo4 = ds.union([cyl1, cyl2])
    bo5 = ds.intersection([cyl1, cyl2])
    # This makes sure the original containers didn't change.
    new_i1 = cyl1["index", "morton_index"]
    new_i1.sort()
    new_i2 = cyl2["index", "morton_index"]
    new_i2.sort()
    assert_array_equal(new_i1, i1)
    assert_array_equal(new_i2, i2)
    # Now make sure the indices also behave as we expect.
    empty = np.array([])
    assert_array_equal(bo1["index", "morton_index"], empty)
    assert_array_equal(bo5["index", "morton_index"], empty)
    b2 = bo2["index", "morton_index"]
    b2.sort()
    assert_array_equal(b2, i1)
    b3 = bo3["index", "morton_index"]
    b3.sort()
    assert_array_equal(b3, ii)
    b4 = bo4["index", "morton_index"]
    b4.sort()
    b5 = bo5["index", "morton_index"]
    b5.sort()
    assert_array_equal(b3, b4)
    bo6 = cyl1 ^ cyl2
    b6 = bo6["index", "morton_index"]
    b6.sort()
    assert_array_equal(b6, np.setxor1d(i1, i2))


def test_boolean_cylinders_overlap():
    r"""Test to make sure that boolean objects (cylinders, overlap)
    behave the way we expect.

    Test overlapping cylinders.
    """
    ds = fake_amr_ds()
    cyl1 = ds.disk([0.45] * 3, [1, 0, 0], 0.2, 0.2)
    cyl2 = ds.disk([0.55] * 3, [1, 0, 0], 0.2, 0.2)
    # Get indices of both.
    i1 = cyl1["index", "morton_index"]
    i2 = cyl2["index", "morton_index"]
    # Make some booleans
    bo1 = cyl1 & cyl2
    bo2 = cyl1 - cyl2
    bo3 = cyl1 | cyl2
    bo4 = ds.union([cyl1, cyl2])
    bo5 = ds.intersection([cyl1, cyl2])
    # Now make sure the indices also behave as we expect.
    vlens = np.intersect1d(i1, i2)
    bite_disk = np.setdiff1d(i1, i2)
    both = np.union1d(i1, i2)
    b1 = bo1["index", "morton_index"]
    b1.sort()
    b2 = bo2["index", "morton_index"]
    b2.sort()
    b3 = bo3["index", "morton_index"]
    b3.sort()
    assert_array_equal(b1, vlens)
    assert_array_equal(b2, bite_disk)
    assert_array_equal(b3, both)
    b4 = bo4["index", "morton_index"]
    b4.sort()
    b5 = bo5["index", "morton_index"]
    b5.sort()
    assert_array_equal(b3, b4)
    assert_array_equal(b1, b5)
    bo6 = cyl1 ^ cyl2
    b6 = bo6["index", "morton_index"]
    b6.sort()
    assert_array_equal(b6, np.setxor1d(i1, i2))
    del ds


def test_boolean_ellipsoids_no_overlap():
    r"""Test to make sure that boolean objects (ellipsoids, no overlap)
    behave the way we expect.

    Test non-overlapping ellipsoids. This also checks that the original
    ellipsoids don't change as part of constructing the booleans.
    """
    ds = fake_amr_ds()
    ell1 = ds.ellipsoid([0.25] * 3, 0.05, 0.05, 0.05, np.array([0.1] * 3), 0.1)
    ell2 = ds.ellipsoid([0.75] * 3, 0.05, 0.05, 0.05, np.array([0.1] * 3), 0.1)
    # Store the original indices
    i1 = ell1["index", "morton_index"]
    i1.sort()
    i2 = ell2["index", "morton_index"]
    i2.sort()
    ii = np.concatenate((i1, i2))
    ii.sort()
    # Make some booleans
    bo1 = ell1 & ell2
    bo2 = ell1 - ell2
    bo3 = ell1 | ell2
    bo4 = ds.union([ell1, ell2])
    bo5 = ds.intersection([ell1, ell2])
    # This makes sure the original containers didn't change.
    new_i1 = ell1["index", "morton_index"]
    new_i1.sort()
    new_i2 = ell2["index", "morton_index"]
    new_i2.sort()
    assert_array_equal(new_i1, i1)
    assert_array_equal(new_i2, i2)
    # Now make sure the indices also behave as we expect.
    empty = np.array([])
    assert_array_equal(bo1["index", "morton_index"], empty)
    assert_array_equal(bo5["index", "morton_index"], empty)
    b2 = bo2["index", "morton_index"]
    b2.sort()
    assert_array_equal(b2, i1)
    b3 = bo3["index", "morton_index"]
    b3.sort()
    assert_array_equal(b3, ii)
    b4 = bo4["index", "morton_index"]
    b4.sort()
    b5 = bo5["index", "morton_index"]
    b5.sort()
    assert_array_equal(b3, b4)
    bo6 = ell1 ^ ell2
    b6 = bo6["index", "morton_index"]
    b6.sort()
    assert_array_equal(b6, np.setxor1d(i1, i2))


def test_boolean_ellipsoids_overlap():
    r"""Test to make sure that boolean objects (ellipsoids, overlap)
    behave the way we expect.

    Test overlapping ellipsoids.
    """
    ds = fake_amr_ds()
    ell1 = ds.ellipsoid([0.45] * 3, 0.05, 0.05, 0.05, np.array([0.1] * 3), 0.1)
    ell2 = ds.ellipsoid([0.55] * 3, 0.05, 0.05, 0.05, np.array([0.1] * 3), 0.1)
    # Get indices of both.
    i1 = ell1["index", "morton_index"]
    i2 = ell2["index", "morton_index"]
    # Make some booleans
    bo1 = ell1 & ell2
    bo2 = ell1 - ell2
    bo3 = ell1 | ell2
    bo4 = ds.union([ell1, ell2])
    bo5 = ds.intersection([ell1, ell2])
    # Now make sure the indices also behave as we expect.
    overlap = np.intersect1d(i1, i2)
    diff = np.setdiff1d(i1, i2)
    both = np.union1d(i1, i2)
    b1 = bo1["index", "morton_index"]
    b1.sort()
    b2 = bo2["index", "morton_index"]
    b2.sort()
    b3 = bo3["index", "morton_index"]
    b3.sort()
    assert_array_equal(b1, overlap)
    assert_array_equal(b2, diff)
    assert_array_equal(b3, both)
    b4 = bo4["index", "morton_index"]
    b4.sort()
    b5 = bo5["index", "morton_index"]
    b5.sort()
    assert_array_equal(b3, b4)
    assert_array_equal(b1, b5)
    bo6 = ell1 ^ ell2
    b6 = bo6["index", "morton_index"]
    b6.sort()
    assert_array_equal(b6, np.setxor1d(i1, i2))


def test_boolean_mix_periodicity():
    r"""Test that a hybrid boolean region behaves as we expect.

    This also tests nested logic and that periodicity works.
    """
    ds = fake_amr_ds()
    re = ds.region([0.5] * 3, [0.0] * 3, [1] * 3)  # whole thing
    sp = ds.sphere([0.95] * 3, 0.3)  # wraps around
    cyl = ds.disk([0.05] * 3, [1, 1, 1], 0.1, 0.4)  # wraps around
    # Get original indices
    rei = re["index", "morton_index"]
    spi = sp["index", "morton_index"]
    cyli = cyl["index", "morton_index"]
    # Make some booleans
    # whole box minux spherical bites at corners
    bo1 = re - sp
    # sphere plus cylinder
    bo2 = sp | cyl
    # a jumble, the region minus the sp+cyl
    bo3 = re - (sp | cyl)
    # Now make sure the indices also behave as we expect.
    bo4 = ds.union([re, sp, cyl])
    bo5 = ds.intersection([re, sp, cyl])
    expect = np.setdiff1d(rei, spi)
    ii = bo1["index", "morton_index"]
    ii.sort()
    assert_array_equal(expect, ii)
    #
    expect = np.union1d(spi, cyli)
    ii = bo2["index", "morton_index"]
    ii.sort()
    assert_array_equal(expect, ii)
    #
    expect = np.union1d(spi, cyli)
    expect = np.setdiff1d(rei, expect)
    ii = bo3["index", "morton_index"]
    ii.sort()
    assert_array_equal(expect, ii)
    b4 = bo4["index", "morton_index"]
    b4.sort()
    b5 = bo5["index", "morton_index"]
    b5.sort()
    ii = np.union1d(np.union1d(rei, cyli), spi)
    ii.sort()
    assert_array_equal(ii, b4)
    ii = np.intersect1d(np.intersect1d(rei, cyli), spi)
    ii.sort()
    assert_array_equal(ii, b5)

    bo6 = (re ^ sp) ^ cyl
    b6 = bo6["index", "morton_index"]
    b6.sort()
    assert_array_equal(b6, np.setxor1d(np.setxor1d(rei, spi), cyli))


def test_boolean_ray_region_no_overlap():
    r"""Test to make sure that boolean objects (ray, region, no overlap)
    behave the way we expect.

    Test non-overlapping ray and region. This also checks that the original
    objects don't change as part of constructing the booleans.
    """
    ds = fake_amr_ds()
    re = ds.box([0.25] * 3, [0.75] * 3)
    ra = ds.ray([0.1] * 3, [0.1, 0.1, 0.9])
    # Store the original indices
    i1 = re["index", "morton_index"]
    i1.sort()
    i2 = ra["index", "morton_index"]
    i2.sort()
    ii = np.concatenate((i1, i2))
    ii.sort()
    # Make some booleans
    bo1 = re & ra
    bo2 = re - ra
    bo3 = re | ra
    bo4 = ds.union([re, ra])
    bo5 = ds.intersection([re, ra])
    # This makes sure the original containers didn't change.
    new_i1 = re["index", "morton_index"]
    new_i1.sort()
    new_i2 = ra["index", "morton_index"]
    new_i2.sort()
    assert_array_equal(new_i1, i1)
    assert_array_equal(new_i2, i2)
    # Now make sure the indices also behave as we expect.
    empty = np.array([])
    assert_array_equal(bo1["index", "morton_index"], empty)
    assert_array_equal(bo5["index", "morton_index"], empty)
    b2 = bo2["index", "morton_index"]
    b2.sort()
    assert_array_equal(b2, i1)
    b3 = bo3["index", "morton_index"]
    b3.sort()
    assert_array_equal(b3, ii)
    b4 = bo4["index", "morton_index"]
    b4.sort()
    b5 = bo5["index", "morton_index"]
    b5.sort()
    assert_array_equal(b3, b4)
    bo6 = re ^ ra
    b6 = bo6["index", "morton_index"]
    b6.sort()
    assert_array_equal(b6, np.setxor1d(i1, i2))


def test_boolean_ray_region_overlap():
    r"""Test to make sure that boolean objects (ray, region, overlap)
    behave the way we expect.

    Test overlapping ray and region. This also checks that the original
    objects don't change as part of constructing the booleans.
    """
    ds = fake_amr_ds()
    re = ds.box([0.25] * 3, [0.75] * 3)
    ra = ds.ray([0] * 3, [1] * 3)
    # Get indices of both.
    i1 = re["index", "morton_index"]
    i2 = ra["index", "morton_index"]
    # Make some booleans
    bo1 = re & ra
    bo2 = re - ra
    bo3 = re | ra
    bo4 = ds.union([re, ra])
    bo5 = ds.intersection([re, ra])
    # Now make sure the indices also behave as we expect.
    short_line = np.intersect1d(i1, i2)
    cube_minus_line = np.setdiff1d(i1, i2)
    both = np.union1d(i1, i2)
    b1 = bo1["index", "morton_index"]
    b1.sort()
    b2 = bo2["index", "morton_index"]
    b2.sort()
    b3 = bo3["index", "morton_index"]
    b3.sort()
    assert_array_equal(b1, short_line)
    assert_array_equal(b2, cube_minus_line)
    assert_array_equal(b3, both)
    b4 = bo4["index", "morton_index"]
    b4.sort()
    b5 = bo5["index", "morton_index"]
    b5.sort()
    assert_array_equal(b3, b4)
    assert_array_equal(b1, b5)
    bo6 = re ^ ra
    b6 = bo6["index", "morton_index"]
    b6.sort()
    assert_array_equal(b6, np.setxor1d(i1, i2))


def test_boolean_rays_no_overlap():
    r"""Test to make sure that boolean objects (rays, no overlap)
    behave the way we expect.

    Test non-overlapping rays.
    """
    ds = fake_amr_ds()
    ra1 = ds.ray([0, 0, 0], [0, 0, 1])
    ra2 = ds.ray([1, 0, 0], [1, 0, 1])
    # Store the original indices
    i1 = ra1["index", "morton_index"]
    i1.sort()
    i2 = ra2["index", "morton_index"]
    i2.sort()
    ii = np.concatenate((i1, i2))
    ii.sort()
    # Make some booleans
    bo1 = ra1 & ra2
    bo2 = ra1 - ra2
    bo3 = ra1 | ra2
    bo4 = ds.union([ra1, ra2])
    bo5 = ds.intersection([ra1, ra2])
    # This makes sure the original containers didn't change.
    new_i1 = ra1["index", "morton_index"]
    new_i1.sort()
    new_i2 = ra2["index", "morton_index"]
    new_i2.sort()
    assert_array_equal(new_i1, i1)
    assert_array_equal(new_i2, i2)
    # Now make sure the indices also behave as we expect.
    empty = np.array([])
    assert_array_equal(bo1["index", "morton_index"], empty)
    assert_array_equal(bo5["index", "morton_index"], empty)
    b2 = bo2["index", "morton_index"]
    b2.sort()
    assert_array_equal(b2, i1)
    b3 = bo3["index", "morton_index"]
    b3.sort()
    assert_array_equal(b3, ii)
    b4 = bo4["index", "morton_index"]
    b4.sort()
    b5 = bo5["index", "morton_index"]
    b5.sort()
    assert_array_equal(b3, b4)
    bo6 = ra1 ^ ra2
    b6 = bo6["index", "morton_index"]
    b6.sort()
    assert_array_equal(b6, np.setxor1d(i1, i2))


def test_boolean_rays_overlap():
    r"""Test to make sure that boolean objects (rays, overlap)
    behave the way we expect.

    Test non-overlapping rays.
    """
    ds = fake_amr_ds()
    ra1 = ds.ray([0] * 3, [1] * 3)
    ra2 = ds.ray([0] * 3, [0.5] * 3)
    # Get indices of both.
    i1 = ra1["index", "morton_index"]
    i1.sort()
    i2 = ra2["index", "morton_index"]
    i2.sort()
    ii = np.concatenate((i1, i2))
    ii.sort()
    # Make some booleans
    bo1 = ra1 & ra2
    bo2 = ra1 - ra2
    bo3 = ra1 | ra2
    bo4 = ds.union([ra1, ra2])
    bo5 = ds.intersection([ra1, ra2])
    # Now make sure the indices also behave as we expect.
    short_line = np.intersect1d(i1, i2)
    short_line_b = np.setdiff1d(i1, i2)
    full_line = np.union1d(i1, i2)
    b1 = bo1["index", "morton_index"]
    b1.sort()
    b2 = bo2["index", "morton_index"]
    b2.sort()
    b3 = bo3["index", "morton_index"]
    b3.sort()
    assert_array_equal(b1, short_line)
    assert_array_equal(b2, short_line_b)
    assert_array_equal(b3, full_line)
    b4 = bo4["index", "morton_index"]
    b4.sort()
    b5 = bo5["index", "morton_index"]
    b5.sort()
    assert_array_equal(b3, i1)
    assert_array_equal(b3, b4)
    assert_array_equal(b1, b5)
    bo6 = ra1 ^ ra2
    b6 = bo6["index", "morton_index"]
    b6.sort()
    assert_array_equal(b6, np.setxor1d(i1, i2))


def test_boolean_slices_no_overlap():
    r"""Test to make sure that boolean objects (slices, no overlap)
    behave the way we expect.

    Test non-overlapping slices. This also checks that the original regions
    don't change as part of constructing the booleans.
    """
    ds = fake_amr_ds()
    sl1 = ds.r[:, :, 0.25]
    sl2 = ds.r[:, :, 0.75]
    # Store the original indices
    i1 = sl1["index", "morton_index"]
    i1.sort()
    i2 = sl2["index", "morton_index"]
    i2.sort()
    ii = np.concatenate((i1, i2))
    ii.sort()
    # Make some booleans
    bo1 = sl1 & sl2
    bo2 = sl1 - sl2
    bo3 = sl1 | sl2
    bo4 = ds.union([sl1, sl2])
    bo5 = ds.intersection([sl1, sl2])
    # This makes sure the original containers didn't change.
    new_i1 = sl1["index", "morton_index"]
    new_i1.sort()
    new_i2 = sl2["index", "morton_index"]
    new_i2.sort()
    assert_array_equal(new_i1, i1)
    assert_array_equal(new_i2, i2)
    # Now make sure the indices also behave as we expect.
    empty = np.array([])
    assert_array_equal(bo1["index", "morton_index"], empty)
    assert_array_equal(bo5["index", "morton_index"], empty)
    b2 = bo2["index", "morton_index"]
    b2.sort()
    assert_array_equal(b2, i1)
    b3 = bo3["index", "morton_index"]
    b3.sort()
    assert_array_equal(b3, ii)
    b4 = bo4["index", "morton_index"]
    b4.sort()
    b5 = bo5["index", "morton_index"]
    b5.sort()
    assert_array_equal(b3, b4)
    bo6 = sl1 ^ sl2
    b6 = bo6["index", "morton_index"]
    b6.sort()
    assert_array_equal(b6, np.setxor1d(i1, i2))


def test_boolean_slices_overlap():
    r"""Test to make sure that boolean objects (slices, overlap)
    behave the way we expect.

    Test overlapping slices.
    """
    ds = fake_amr_ds()
    sl1 = ds.r[:, :, 0.25]
    sl2 = ds.r[:, 0.75, :]
    # Get indices of both.
    i1 = sl1["index", "morton_index"]
    i2 = sl2["index", "morton_index"]
    # Make some booleans
    bo1 = sl1 & sl2
    bo2 = sl1 - sl2
    bo3 = sl1 | sl2
    bo4 = ds.union([sl1, sl2])
    bo5 = ds.intersection([sl1, sl2])
    # Now make sure the indices also behave as we expect.
    line = np.intersect1d(i1, i2)
    orig = np.setdiff1d(i1, i2)
    both = np.union1d(i1, i2)
    b1 = bo1["index", "morton_index"]
    b1.sort()
    b2 = bo2["index", "morton_index"]
    b2.sort()
    b3 = bo3["index", "morton_index"]
    b3.sort()
    assert_array_equal(b1, line)
    assert_array_equal(b2, orig)
    assert_array_equal(b3, both)
    b4 = bo4["index", "morton_index"]
    b4.sort()
    b5 = bo5["index", "morton_index"]
    b5.sort()
    assert_array_equal(b3, b4)
    assert_array_equal(b1, b5)
    bo6 = sl1 ^ sl2
    b6 = bo6["index", "morton_index"]
    b6.sort()
    assert_array_equal(b6, np.setxor1d(i1, i2))


def test_boolean_ray_slice_no_overlap():
    r"""Test to make sure that boolean objects (ray, slice, no overlap)
    behave the way we expect.

    Test non-overlapping ray and slice. This also checks that the original
    regions don't change as part of constructing the booleans.
    """
    ds = fake_amr_ds()
    sl = ds.r[:, :, 0.25]
    ra = ds.ray([0] * 3, [0, 1, 0])
    # Store the original indices
    i1 = sl["index", "morton_index"]
    i1.sort()
    i2 = ra["index", "morton_index"]
    i2.sort()
    ii = np.concatenate((i1, i2))
    ii.sort()
    # Make some booleans
    bo1 = sl & ra
    bo2 = sl - ra
    bo3 = sl | ra
    bo4 = ds.union([sl, ra])
    bo5 = ds.intersection([sl, ra])
    # This makes sure the original containers didn't change.
    new_i1 = sl["index", "morton_index"]
    new_i1.sort()
    new_i2 = ra["index", "morton_index"]
    new_i2.sort()
    assert_array_equal(new_i1, i1)
    assert_array_equal(new_i2, i2)
    # Now make sure the indices also behave as we expect.
    empty = np.array([])
    assert_array_equal(bo1["index", "morton_index"], empty)
    assert_array_equal(bo5["index", "morton_index"], empty)
    b2 = bo2["index", "morton_index"]
    b2.sort()
    assert_array_equal(b2, i1)
    b3 = bo3["index", "morton_index"]
    b3.sort()
    assert_array_equal(b3, ii)
    b4 = bo4["index", "morton_index"]
    b4.sort()
    b5 = bo5["index", "morton_index"]
    b5.sort()
    assert_array_equal(b3, b4)
    bo6 = sl ^ ra
    b6 = bo6["index", "morton_index"]
    b6.sort()
    assert_array_equal(b6, np.setxor1d(i1, i2))


def test_boolean_ray_slice_overlap():
    r"""Test to make sure that boolean objects (rays and slices, overlap)
    behave the way we expect.

    Test overlapping rays and slices.
    """
    ds = fake_amr_ds()
    sl = ds.r[:, :, 0.25]
    ra = ds.ray([0, 0, 0.25], [0, 1, 0.25])
    # Get indices of both.
    i1 = sl["index", "morton_index"]
    i1.sort()
    i2 = ra["index", "morton_index"]
    i1.sort()
    ii = np.concatenate((i1, i2))
    ii.sort()
    # Make some booleans
    bo1 = sl & ra
    bo2 = sl - ra
    bo3 = sl | ra
    bo4 = ds.union([sl, ra])
    bo5 = ds.intersection([sl, ra])
    # Now make sure the indices also behave as we expect.
    line = np.intersect1d(i1, i2)
    sheet_minus_line = np.setdiff1d(i1, i2)
    sheet = np.union1d(i1, i2)
    b1 = bo1["index", "morton_index"]
    b1.sort()
    b2 = bo2["index", "morton_index"]
    b2.sort()
    b3 = bo3["index", "morton_index"]
    b3.sort()
    assert_array_equal(b1, line)
    assert_array_equal(b2, sheet_minus_line)
    assert_array_equal(b3, sheet)
    b4 = bo4["index", "morton_index"]
    b4.sort()
    b5 = bo5["index", "morton_index"]
    b5.sort()
    assert_array_equal(b3, i1)
    assert_array_equal(b3, b4)
    assert_array_equal(b1, b5)
    bo6 = sl ^ ra
    b6 = bo6["index", "morton_index"]
    b6.sort()
    assert_array_equal(b6, np.setxor1d(i1, i2))
