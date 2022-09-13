import numpy as np

from yt.testing import assert_almost_equal, assert_equal, fake_sph_grid_ds

n_ref = 4


def test_building_tree():
    """
    Build an octree and make sure correct number of particles
    """
    ds = fake_sph_grid_ds()
    octree = ds.octree(n_ref=n_ref)
    assert octree[("index", "x")].shape[0] == 17


def test_sph_interpolation_scatter():
    """
    Generate an octree, perform some SPH interpolation and check with some
    answer testing
    """

    ds = fake_sph_grid_ds(hsml_factor=26.0)
    ds._sph_ptypes = ("io",)

    ds.use_sph_normalization = False

    octree = ds.octree(n_ref=n_ref)

    density = octree[("io", "density")]
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

    density = octree[("io", "density")]
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

    depth = octree[("index", "depth")]
    depth_ans = np.array([0] + [1] * 8 + [2] * 8, dtype=np.int64)
    assert_equal(depth, depth_ans)

    size_ans = np.zeros((depth.shape[0], 3), dtype=np.float64)
    for i in range(size_ans.shape[0]):
        size_ans[i, :] = (ds.domain_right_edge - ds.domain_left_edge) / 2.0 ** depth[i]

    dx = octree[("index", "dx")].d
    assert_almost_equal(dx, size_ans[:, 0])

    dy = octree[("index", "dy")].d
    assert_almost_equal(dy, size_ans[:, 1])

    dz = octree[("index", "dz")].d
    assert_almost_equal(dz, size_ans[:, 2])

    refined = octree[("index", "refined")]
    refined_ans = np.array([True] + [False] * 7 + [True] + [False] * 8, dtype=np.bool_)
    assert_equal(refined, refined_ans)
