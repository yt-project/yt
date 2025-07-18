from typing import Any, Optional

import numpy as np

from yt.data_objects.selection_objects.region import YTRegion
from yt.geometry.selection_routines import (
    bbox_intersects,
    fully_contains,
)
from yt.utilities.lib.geometry_utils import get_hilbert_indices

# State diagram to compute the hilbert curve
_STATE_DIAGRAM = np.array(
    [
        [
            [1, 2, 0, 6, 11, 4, 5, 6, 10, 4, 7, 10],
            [0, 0, 0, 2, 4, 6, 4, 6, 2, 2, 4, 6],
        ],
        [
            [2, 6, 9, 0, 11, 4, 7, 1, 3, 4, 2, 3],
            [1, 7, 3, 3, 3, 5, 7, 7, 5, 1, 5, 1],
        ],
        [
            [3, 0, 10, 6, 0, 8, 5, 6, 1, 8, 11, 2],
            [3, 1, 7, 1, 5, 1, 3, 5, 3, 5, 7, 7],
        ],
        [
            [2, 7, 9, 11, 7, 8, 3, 10, 1, 8, 2, 6],
            [2, 6, 4, 0, 2, 2, 0, 4, 4, 6, 6, 0],
        ],
        [
            [4, 8, 1, 9, 5, 0, 1, 9, 10, 2, 7, 10],
            [7, 3, 1, 5, 7, 7, 5, 1, 1, 3, 3, 5],
        ],
        [
            [5, 8, 1, 0, 9, 6, 1, 4, 3, 7, 5, 3],
            [6, 4, 2, 4, 0, 4, 6, 0, 6, 0, 2, 2],
        ],
        [
            [3, 0, 11, 9, 0, 10, 11, 9, 5, 2, 8, 4],
            [4, 2, 6, 6, 6, 0, 2, 2, 0, 4, 0, 4],
        ],
        [
            [5, 7, 11, 8, 7, 6, 11, 10, 9, 3, 5, 4],
            [5, 5, 5, 7, 1, 3, 1, 3, 7, 7, 1, 3],
        ],
    ]
)


def hilbert3d(
    ijk: "np.ndarray[Any, np.dtype[np.int64]]", bit_length: int
) -> "np.ndarray[Any, np.dtype[np.float64]]":
    """Compute the order using Hilbert indexing.

    Arguments
    ---------
    ijk : (N, ndim) integer array
      The positions
    bit_length : integer
      The bit_length for the indexing.
    """
    ijk = np.atleast_2d(ijk)
    # A note here: there is a freedom in the way hilbert indices are
    # being computed (should it be xyz or yzx or zxy etc.)
    # and the yt convention is not the same as the RAMSES one.
    return get_hilbert_indices(bit_length, ijk[:, [1, 2, 0]].astype(np.int64))


def get_intersecting_cpus(
    ds,
    region: YTRegion,
    LE: Optional["np.ndarray[Any, np.dtype[np.float64]]"] = None,
    dx: float = 1.0,
    dx_cond: float | None = None,
    factor: float = 4.0,
    bound_keys: Optional["np.ndarray[Any, np.dtype[np.float64]]"] = None,
) -> set[int]:
    """
    Find the subset of CPUs that intersect the bbox in a recursive fashion.
    """
    if LE is None:
        LE = np.array([0, 0, 0], dtype="d")
    if dx_cond is None:
        bbox = region.get_bbox()
        dx_cond = float((bbox[1] - bbox[0]).min().to("code_length"))
    if bound_keys is None:
        ncpu = ds.parameters["ncpu"]
        bound_keys = np.empty(ncpu + 1, dtype="float64")
        bound_keys[:ncpu] = [ds.hilbert_indices[icpu + 1][0] for icpu in range(ncpu)]
        bound_keys[ncpu] = ds.hilbert_indices[ncpu][1]

    # If the current dx is smaller than the smallest size of the bbox
    if dx < dx_cond / factor:
        # Finish recursion
        return get_cpu_list_cuboid(ds, np.asarray([LE, LE + dx]), bound_keys)

    # If the current cell is fully within the selected region, stop recursion
    if fully_contains(region.selector, LE, dx):
        return get_cpu_list_cuboid(ds, np.asarray([LE, LE + dx]), bound_keys)

    dx /= 2

    ret = set()
    # Compute intersection of the eight subcubes with the bbox and recurse.
    for i in range(2):
        for j in range(2):
            for k in range(2):
                LE_new = LE + np.array([i, j, k], dtype="d") * dx

                if bbox_intersects(region.selector, LE_new, dx):
                    ret.update(
                        get_intersecting_cpus(
                            ds, region, LE_new, dx, dx_cond, factor, bound_keys
                        )
                    )
    return ret


def get_cpu_list_cuboid(
    ds,
    X: "np.ndarray[Any, np.dtype[np.float64]]",
    bound_keys: "np.ndarray[Any, np.dtype[np.float64]]",
) -> set[int]:
    """
    Return the list of the CPU intersecting with the cuboid containing the positions.
    Note that it will be 0-indexed.

    Parameters
    ----------
    ds : Dataset
      The dataset containing the information
    X : (N, ndim) float array
      An array containing positions. They should be between 0 and 1.
    """
    X = np.atleast_2d(X)
    if X.shape[1] != 3:
        raise NotImplementedError("This function is only implemented in 3D.")

    levelmax = ds.parameters["levelmax"]
    ndim = ds.parameters["ndim"]

    xmin, ymin, zmin = X.min(axis=0)
    xmax, ymax, zmax = X.max(axis=0)

    dmax = max(xmax - xmin, ymax - ymin, zmax - zmin)
    ilevel = int(np.ceil(-np.log2(dmax)))

    lmin = ilevel
    bit_length = lmin - 1
    maxdom = 2**bit_length

    imin, imax, jmin, jmax, kmin, kmax = 0, 0, 0, 0, 0, 0
    if bit_length > 0:
        imin = int(xmin * maxdom)
        imax = imin + 1
        jmin = int(ymin * maxdom)
        jmax = jmin + 1
        kmin = int(zmin * maxdom)
        kmax = kmin + 1

    dkey = (2 ** (levelmax + 1) / maxdom) ** ndim
    ndom = 1
    if bit_length > 0:
        ndom = 8

    ijkdom = idom, jdom, kdom = np.empty((3, 8), dtype=np.int64)

    idom[0], idom[1] = imin, imax
    idom[2], idom[3] = imin, imax
    idom[4], idom[5] = imin, imax
    idom[6], idom[7] = imin, imax

    jdom[0], jdom[1] = jmin, jmin
    jdom[2], jdom[3] = jmax, jmax
    jdom[4], jdom[5] = jmin, jmin
    jdom[6], jdom[7] = jmax, jmax

    kdom[0], kdom[1] = kmin, kmin
    kdom[2], kdom[3] = kmin, kmin
    kdom[4], kdom[5] = kmax, kmax
    kdom[6], kdom[7] = kmax, kmax

    bounding_min, bounding_max = np.zeros(ndom), np.zeros(ndom)
    if bit_length > 0:
        order_min = hilbert3d(ijkdom.T, bit_length)
    for i in range(ndom):
        if bit_length > 0:
            omin = order_min[i]
        else:
            omin = 0
        bounding_min[i] = omin * dkey
        bounding_max[i] = (omin + 1) * dkey

    cpu_min = np.searchsorted(bound_keys, bounding_min, side="right") - 1
    cpu_max = np.searchsorted(bound_keys, bounding_max, side="right")

    cpu_read: set[int] = set()

    for i in range(ndom):
        cpu_read.update(range(cpu_min[i], cpu_max[i]))

    return cpu_read
