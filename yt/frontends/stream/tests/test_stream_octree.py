import yt
import numpy as np

OCT_MASK_LIST = [8, 0, 0, 0, 0, 8, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0,
                 8, 0, 0, 0, 0, 0, 0, 0,
                 0]


def test_octree():
    # See Issue #1272
    octree_mask = np.array(OCT_MASK_LIST, dtype=np.uint8)

    quantities = {}
    quantities[('gas', 'density')] = np.ones((22, 1), dtype=float)

    bbox = np.array([[-10., 10.], [-10., 10.], [-10., 10.]])

    ds = yt.load_octree(octree_mask=octree_mask,
                        data=quantities,
                        bbox=bbox,
                        over_refine_factor=0,
                        partial_coverage=0)

    proj = ds.proj('density', 'x')
    proj['density']
