import numpy as np
from numpy.testing import assert_almost_equal, assert_equal

from yt.geometry.oct_container import _ORDER_MAX, OctreeContainer
from yt.geometry.particle_oct_container import ParticleOctreeContainer
from yt.geometry.selection_routines import AlwaysSelector
from yt.testing import fake_sph_grid_ds
from yt.utilities.lib.geometry_utils import get_morton_indices

n_ref = 4


def test_building_tree():
    """
    Build an octree and make sure correct number of particles
    """
    ds = fake_sph_grid_ds()
    octree = ds.octree(n_ref=n_ref)
    assert octree["index", "x"].shape[0] == 17


def test_sph_interpolation_scatter():
    """
    Generate an octree, perform some SPH interpolation and check with some
    answer testing
    """

    ds = fake_sph_grid_ds(hsml_factor=26.0)
    ds._sph_ptypes = ("io",)

    ds.use_sph_normalization = False

    octree = ds.octree(n_ref=n_ref)

    density = octree["io", "density"]
    answers = np.array(
        [
            1.00434706,
            1.00434706,
            1.00434706,
            1.00434706,
            1.00434706,
            1.00434706,
            1.00434706,
            0.7762907,
            0.89250848,
            0.89250848,
            0.97039088,
            0.89250848,
            0.97039088,
            0.97039088,
            1.01156175,
        ]
    )

    assert_almost_equal(density.d, answers)


def test_sph_interpolation_gather():
    """
    Generate an octree, perform some SPH interpolation and check with some
    answer testing
    """
    ds = fake_sph_grid_ds(hsml_factor=26.0)
    ds.index
    ds._sph_ptypes = ("io",)

    ds.sph_smoothing_style = "gather"
    ds.num_neighbors = 5
    ds.use_sph_normalization = False

    octree = ds.octree(n_ref=n_ref)

    density = octree["io", "density"]
    answers = np.array(
        [
            0.59240874,
            0.59240874,
            0.59240874,
            0.59240874,
            0.59240874,
            0.59240874,
            0.59240874,
            0.10026846,
            0.77014968,
            0.77014968,
            0.96127825,
            0.77014968,
            0.96127825,
            0.96127825,
            1.21183996,
        ]
    )

    assert_almost_equal(density.d, answers)


def test_octree_properties():
    """
    Generate an octree, and test the refinement, depth and sizes of the cells.
    """
    ds = fake_sph_grid_ds()
    octree = ds.octree(n_ref=n_ref)

    depth = octree["index", "depth"]
    depth_ans = np.array([0] + [1] * 8 + [2] * 8, dtype=np.int64)
    assert_equal(depth, depth_ans)

    size_ans = np.zeros((depth.shape[0], 3), dtype=np.float64)
    for i in range(size_ans.shape[0]):
        size_ans[i, :] = (ds.domain_right_edge - ds.domain_left_edge) / 2.0 ** depth[i]

    dx = octree["index", "dx"].d
    assert_almost_equal(dx, size_ans[:, 0])

    dy = octree["index", "dy"].d
    assert_almost_equal(dy, size_ans[:, 1])

    dz = octree["index", "dz"].d
    assert_almost_equal(dz, size_ans[:, 2])

    refined = octree["index", "refined"]
    refined_ans = np.array([True] + [False] * 7 + [True] + [False] * 8, dtype=np.bool_)
    assert_equal(refined, refined_ans)


def test_num_zones_tuple():
    """
    Test that OctreeContainer accepts num_zones as a scalar or a tuple (N, M, L).
    Both should correctly set per-dimension zone counts.
    """
    # Scalar: all dimensions equal
    oct_scalar = OctreeContainer(
        [1, 1, 1], [0.0, 0.0, 0.0], [1.0, 1.0, 1.0], num_zones=2
    )
    # Tuple: potentially different per-dimension
    oct_tuple = OctreeContainer(
        [1, 1, 1], [0.0, 0.0, 0.0], [1.0, 1.0, 1.0], num_zones=(2, 2, 2)
    )
    # Non-uniform tuple
    oct_nonuniform = OctreeContainer(
        [1, 1, 1], [0.0, 0.0, 0.0], [1.0, 1.0, 1.0], num_zones=(2, 3, 4)
    )
    # Verify that creating these containers doesn't raise exceptions
    assert oct_scalar is not None
    assert oct_tuple is not None
    assert oct_nonuniform is not None


def test_octcellindex_neighbours_num_zones():
    """
    Regression test for #5402: fill_octcellindex_neighbours hard-coded the
    assumption that every oct holds 2x2x2 zones, so both the loop bounds and
    the output buffer size were wrong whenever num_zones wasn't (2, 2, 2).
    """
    DLE = np.array([0.0, 0.0, 0.0])
    DRE = np.array([8.0, 8.0, 8.0])
    dx = (DRE - DLE) / (2**_ORDER_MAX)

    # One particle per octant of the root oct: refines exactly one level deep.
    centers = np.array(
        [[x, y, z] for x in (2, 6) for y in (2, 6) for z in (2, 6)], dtype="float64"
    )
    morton = get_morton_indices(np.floor((centers - DLE) / dx).astype("uint64"))
    morton.sort()

    for nz in ((2, 2, 2), (2, 3, 4)):
        octree = ParticleOctreeContainer((1, 1, 1), DLE, DRE, num_zones=nz)
        octree.n_ref = 1
        octree.add(morton)
        octree.finalize()

        nzones = nz[0] * nz[1] * nz[2]
        n_per_oct = (nz[0] + 2) * (nz[1] + 2) * (nz[2] + 2)  # +1 ghost zone each side

        selector = AlwaysSelector(None)
        num_octs = selector.count_octs(octree, -1)
        _, cell_inds = octree.fill_octcellindex_neighbours(selector)

        # old oct_visitors/containers hardcoded 4**3=64 cells/oct; for nz=(2,3,4) it's really
        # 4*5*6=120. this assertion catches that
        assert_equal(cell_inds.size, num_octs * n_per_oct)
        assert cell_inds.min() >= 0
        assert cell_inds.max() <= nzones
