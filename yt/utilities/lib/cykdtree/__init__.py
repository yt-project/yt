from yt.utilities.lib.cykdtree import plot  # NOQA
from yt.utilities.lib.cykdtree.kdtree import PyKDTree, PyNode  # NOQA


def make_tree(pts, **kwargs):
    r"""Build a KD-tree for a set of points.

    Args:
        pts (np.ndarray of float64): (n,m) Array of n mD points.
        \*\*kwargs: Additional keyword arguments are passed to the appropriate
            class for constructuing the tree.

    Returns:
        T (:class:`cykdtree.PyKDTree`): KDTree object.

    Raises:
        ValueError: If `pts` is not a 2D array.

    """
    # Check input
    if pts.ndim != 2:
        raise ValueError("pts must be a 2D array of ND coordinates")
    T = PyKDTree(pts, **kwargs)
    return T
