from yt.testing import *
from yt.data_objects.api import add_field

def setup():
    from yt.config import ytcfg
    ytcfg["yt","__withintesting"] = "True"
    def _ID(field, data):
        width = data.pf.domain_right_edge - data.pf.domain_left_edge
        min_dx = 1.0/8192
        delta = width / min_dx
        x = data['x'] - min_dx / 2.
        y = data['y'] - min_dx / 2.
        z = data['z'] - min_dx / 2.
        xi = x / min_dx
        yi = y / min_dx
        zi = z / min_dx
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
    for n in [1, 2, 4, 8]:
        pf = fake_random_pf(64, nprocs=n)
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
        yield assert_array_equal, new_i1, i1
        yield assert_array_equal, new_i2, i2
        # Now make sure the indices also behave as we expect.
        empty = np.array([])
        yield assert_array_equal, bo1['ID'], empty
        b2 = bo2['ID']
        b2.sort()
        yield assert_array_equal, b2, i1
        b3 = bo3['ID']
        b3.sort()
        yield assert_array_equal, b3, ii
 
def test_boolean_spheres_overlap():
    r"""Test to make sure that boolean objects (spheres, overlap)
    behave the way we expect.

    Test overlapping spheres.
    """
    for n in [1, 2, 4, 8]:
        pf = fake_random_pf(64, nprocs=n)
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
        yield assert_array_equal, b1, lens
        yield assert_array_equal, b2, apple
        yield assert_array_equal, b3, both

def test_boolean_regions_no_overlap():
    r"""Test to make sure that boolean objects (regions, no overlap)
    behave the way we expect.

    Test non-overlapping regions. This also checks that the original regions
    don't change as part of constructing the booleans.
    """
    for n in [1, 2, 4, 8]:
        pf = fake_random_pf(64, nprocs=n)
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
        yield assert_array_equal, new_i1, i1
        yield assert_array_equal, new_i2, i2
        # Now make sure the indices also behave as we expect.
        empty = np.array([])
        yield assert_array_equal, bo1['ID'], empty
        b2 = bo2['ID']
        b2.sort()
        yield assert_array_equal, b2, i1 
        b3 = bo3['ID']
        b3.sort()
        yield assert_array_equal, b3, ii

def test_boolean_regions_overlap():
    r"""Test to make sure that boolean objects (regions, overlap)
    behave the way we expect.

    Test overlapping regions.
    """
    for n in [1, 2, 4, 8]:
        pf = fake_random_pf(64, nprocs=n)
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
        yield assert_array_equal, b1, cube
        yield assert_array_equal, b2, bite_cube
        yield assert_array_equal, b3, both

def test_boolean_cylinders_no_overlap():
    r"""Test to make sure that boolean objects (cylinders, no overlap)
    behave the way we expect.

    Test non-overlapping cylinders. This also checks that the original cylinders
    don't change as part of constructing the booleans.
    """
    for n in [1, 2, 4, 8]:
        pf = fake_random_pf(64, nprocs=n)
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
        yield assert_array_equal, new_i1, i1
        yield assert_array_equal, new_i2, i2
        # Now make sure the indices also behave as we expect.
        empty = np.array([])
        yield assert_array_equal, bo1['ID'], empty
        b2 = bo2['ID']
        b2.sort()
        yield assert_array_equal, b2, i1
        b3 = bo3['ID']
        b3.sort()
        yield assert_array_equal, b3, ii

def test_boolean_cylinders_overlap():
    r"""Test to make sure that boolean objects (cylinders, overlap)
    behave the way we expect.

    Test overlapping cylinders.
    """
    for n in [1, 2, 4, 8]:
        pf = fake_random_pf(64, nprocs=n)
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
        yield assert_array_equal, b1, vlens
        yield assert_array_equal, b2, bite_disk
        yield assert_array_equal, b3, both

def test_boolean_ellipsoids_no_overlap():
    r"""Test to make sure that boolean objects (ellipsoids, no overlap)
    behave the way we expect.

    Test non-overlapping ellipsoids. This also checks that the original
    ellipsoids don't change as part of constructing the booleans.
    """
    for n in [1, 2, 4, 8]:
        pf = fake_random_pf(64, nprocs=n)
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
        yield assert_array_equal, new_i1, i1 
        yield assert_array_equal, new_i2, i2
        # Now make sure the indices also behave as we expect.
        empty = np.array([])
        yield assert_array_equal, bo1['ID'], empty
        b2 = bo2['ID']
        b2.sort()
        yield assert_array_equal, b2, i1
        b3 = bo3['ID']
        b3.sort()
        yield assert_array_equal, b3, ii

def test_boolean_ellipsoids_overlap():
    r"""Test to make sure that boolean objects (ellipsoids, overlap)
    behave the way we expect.

    Test overlapping ellipsoids.
    """
    for n in [1, 2, 4, 8]:
        pf = fake_random_pf(64, nprocs=n)
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
        yield assert_array_equal, b1, overlap
        yield assert_array_equal, b2, diff
        yield assert_array_equal, b3, both

def test_boolean_mix_periodicity():
    r"""Test that a hybrid boolean region behaves as we expect.

    This also tests nested logic and that periodicity works.
    """
    for n in [1, 2, 4, 8]:
        pf = fake_random_pf(64, nprocs=n)
        pf.h
        re = pf.h.region([0.5]*3, [0.0]*3, [1]*3) # whole thing
        sp = pf.h.sphere([0.95]*3, 0.3) # wraps around
        cyl = pf.h.disk([0.05]*3, [1,1,1], 0.1, 0.4) # wraps around
        # Get original indices
        rei = re['ID']
        spi = sp['ID']
        cyli = cyl['ID']
        # Make some booleans
        # whole box minux spherical bites at corners
        bo1 = pf.h.boolean([re, "NOT", sp])
        # sphere plus cylinder
        bo2 = pf.h.boolean([sp, "OR", cyl])
        # a jumble, the region minus the sp+cyl
        bo3 = pf.h.boolean([re, "NOT", "(", sp, "OR", cyl, ")"])
        # Now make sure the indices also behave as we expect.
        expect = np.setdiff1d(rei, spi)
        ii = bo1['ID']
        ii.sort()
        yield assert_array_equal, expect, ii
        #
        expect = np.union1d(spi, cyli)
        ii = bo2['ID']
        ii.sort()
        yield assert_array_equal, expect, ii
        #
        expect = np.union1d(spi, cyli)
        expect = np.setdiff1d(rei, expect)
        ii = bo3['ID']
        ii.sort()
        yield assert_array_equal, expect, ii

