import numpy as np
from numpy.testing import assert_equal

import yt

OCT_MASK_LIST = [
    8,
    0,
    0,
    0,
    0,
    8,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    8,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
]


def test_octree():
    # See Issue #1272
    octree_mask = np.array(OCT_MASK_LIST, dtype=np.uint8)

    quantities = {}
    quantities[("gas", "density")] = np.random.random((22, 1))

    bbox = np.array([[-10.0, 10.0], [-10.0, 10.0], [-10.0, 10.0]])

    ds = yt.load_octree(
        octree_mask=octree_mask,
        data=quantities,
        bbox=bbox,
        num_zones=1,
        partial_coverage=0,
    )

    proj = ds.proj(("gas", "density"), "x")
    proj[("gas", "density")]

    assert_equal(ds.r[:]["ones"].size, 22)
    rho1 = quantities["gas", "density"].ravel()
    rho2 = ds.r[:]["density"].copy()
    rho1.sort()
    rho2.sort()
    assert_equal(rho1, rho2)
