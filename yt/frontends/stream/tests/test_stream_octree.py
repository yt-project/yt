import numpy as np

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
    quantities[("gas", "density")] = np.ones((22, 1), dtype="float64")

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
