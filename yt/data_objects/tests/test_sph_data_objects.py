import numpy as np

from yt.testing import \
    assert_equal, \
    fake_sph_orientation_ds

# The number of particles along each slice axis at that coordinate
SLICE_ANSWERS = {
    ('x', 0): 6,
    ('x', 0.5): 0,
    ('x', 1): 1,
    ('y', 0): 5,
    ('y', 1): 1,
    ('y', 2): 1,
    ('z', 0): 4,
    ('z', 1): 1,
    ('z', 2): 1,
    ('z', 3): 1,
}

def test_slice():
    ds = fake_sph_orientation_ds()
    for (ax, coord), answer in SLICE_ANSWERS.items():
        # test that we can still select particles even if we offset the slice
        # within each particle's smoothing volumes
        for i in range(-1, 2):
            sl = ds.slice(ax, coord + i*0.1)
            assert_equal(sl['gas', 'density'].shape[0], answer)

REGION_ANSWERS = {
    ((-4, -4, -4), (4, 4, 4)): 7,
    ((0, 0, 0), (4, 4, 4)): 7,
    ((1, 0, 0), (4, 4, 4)): 1,
    ((0, 1, 0), (4, 4, 4)): 2,
    ((0, 0, 1), (4, 4, 4)): 3,
    ((0, 0, 0), (4, 4, 2)): 6,
    ((0, 0, 0), (4, 4, 1)): 5,
    ((0, 0, 0), (4, 1, 4)): 6,
    ((0, 0, 0), (1, 1, 4)): 6,
}

def test_region():
    ds = fake_sph_orientation_ds()
    for (left_edge, right_edge), answer in REGION_ANSWERS.items():
        # test that regions enclosing a particle's smoothing region
        # correctly select SPH particles
        for i in range(-1, 2):
            left_edge = np.array([le + i*0.1 for le in left_edge])
            right_edge = np.array([re + i*0.1 for re in right_edge])
            # check if we went off the edge of the domain
            whl = left_edge < ds.domain_left_edge
            left_edge[whl] = ds.domain_left_edge[whl]
            whr = left_edge < ds.domain_left_edge
            left_edge[whr] = ds.domain_left_edge[whr]
            reg = ds.box(left_edge, right_edge)
            assert_equal(reg['gas', 'density'].shape[0], answer)

