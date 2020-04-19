import numpy as np

from yt.data_objects.construction_data_containers import YTOctree
from yt.testing import assert_almost_equal, fake_sph_grid_ds

n_ref = 4

def test_building_tree():
    """
    Test function to build an octree and make sure correct number of particles
    """
    ds = fake_sph_grid_ds()
    octree = ds.octree(n_ref=n_ref)
    assert(type(octree) == YTOctree)
    assert(octree[('index', 'x')].shape[0] == 17)


def test_sph_interpolation_scatter():
    """
    Just generate an octree, perform some SPH interpolation and check with some
    answer testing
    """

    ds = fake_sph_grid_ds(hsml_factor=26.0)
    ds._sph_ptypes = ('io',)

    octree = ds.octree(n_ref=n_ref)
    octree.use_sph_normalization = False

    density = octree[('io', 'density')]
    answers = np.array([
        1.00434706, 1.00434706, 1.00434706, 1.00434706,
        1.00434706, 1.00434706, 1.00434706, 0.7762907,
        0.89250848, 0.89250848, 0.97039088, 0.89250848,
        0.97039088, 0.97039088, 1.01156175])

    assert_almost_equal(density.d, answers)


def test_sph_interpolation_gather():
    """
    Just generate an octree, perform some SPH interpolation and check with some
    answer testing
    """
    ds = fake_sph_grid_ds(hsml_factor=26.0)
    ds.index
    ds._sph_ptypes = ('io',)

    octree = ds.octree(n_ref=n_ref)
    octree.sph_smoothing_style = 'gather'
    octree.num_neighbors = 5
    octree.use_sph_normalization = False
    print(octree._sph_smoothing_style)

    density = octree[('io', 'density')]
    answers = np.array([0.59240874, 0.59240874, 0.59240874, 0.59240874,
                        0.59240874, 0.59240874, 0.59240874, 0.10026846,
                        0.77014968, 0.77014968, 0.96127825, 0.77014968,
                        0.96127825, 0.96127825, 1.21183996])

    assert_almost_equal(density.d, answers)
