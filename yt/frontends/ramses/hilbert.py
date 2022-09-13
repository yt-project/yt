import numpy as np


def hilbert3d(X, bit_length):
    """Compute the order using Hilbert indexing.

    Arguments
    ---------
    X : (N, ndim) float array
      The positions
    bit_length : integer
      The bit_length for the indexing.
    """
    X = np.atleast_2d(X)
    state_diagram = (
        np.array(
            [
                1,
                2,
                3,
                2,
                4,
                5,
                3,
                5,
                0,
                1,
                3,
                2,
                7,
                6,
                4,
                5,
                2,
                6,
                0,
                7,
                8,
                8,
                0,
                7,
                0,
                7,
                1,
                6,
                3,
                4,
                2,
                5,
                0,
                9,
                10,
                9,
                1,
                1,
                11,
                11,
                0,
                3,
                7,
                4,
                1,
                2,
                6,
                5,
                6,
                0,
                6,
                11,
                9,
                0,
                9,
                8,
                2,
                3,
                1,
                0,
                5,
                4,
                6,
                7,
                11,
                11,
                0,
                7,
                5,
                9,
                0,
                7,
                4,
                3,
                5,
                2,
                7,
                0,
                6,
                1,
                4,
                4,
                8,
                8,
                0,
                6,
                10,
                6,
                6,
                5,
                1,
                2,
                7,
                4,
                0,
                3,
                5,
                7,
                5,
                3,
                1,
                1,
                11,
                11,
                4,
                7,
                3,
                0,
                5,
                6,
                2,
                1,
                6,
                1,
                6,
                10,
                9,
                4,
                9,
                10,
                6,
                7,
                5,
                4,
                1,
                0,
                2,
                3,
                10,
                3,
                1,
                1,
                10,
                3,
                5,
                9,
                2,
                5,
                3,
                4,
                1,
                6,
                0,
                7,
                4,
                4,
                8,
                8,
                2,
                7,
                2,
                3,
                2,
                1,
                5,
                6,
                3,
                0,
                4,
                7,
                7,
                2,
                11,
                2,
                7,
                5,
                8,
                5,
                4,
                5,
                7,
                6,
                3,
                2,
                0,
                1,
                10,
                3,
                2,
                6,
                10,
                3,
                4,
                4,
                6,
                1,
                7,
                0,
                5,
                2,
                4,
                3,
            ]
        )
        .reshape(12, 2, 8)
        .T
    )

    x_bit_mask, y_bit_mask, z_bit_mask = (
        np.zeros(bit_length, dtype=bool) for _ in range(3)
    )
    i_bit_mask = np.zeros(3 * bit_length, dtype=bool)

    npoint = X.shape[0]
    order = np.zeros(npoint)

    # Convert positions to binary
    for ip in range(npoint):
        for i in range(bit_length):
            mask = 0b01 << i
            x_bit_mask[i] = X[ip, 0] & mask
            y_bit_mask[i] = X[ip, 1] & mask
            z_bit_mask[i] = X[ip, 2] & mask

        for i in range(bit_length):
            # Interleave bits
            i_bit_mask[3 * i + 2] = x_bit_mask[i]
            i_bit_mask[3 * i + 1] = y_bit_mask[i]
            i_bit_mask[3 * i] = z_bit_mask[i]

        # Build Hilbert ordering using state diagram
        cstate = 0
        for i in range(bit_length - 1, -1, -1):
            sdigit = (
                4 * i_bit_mask[3 * i + 2]
                + 2 * i_bit_mask[3 * i + 1]
                + 1 * i_bit_mask[3 * i]
            )
            nstate = state_diagram[sdigit, 0, cstate]
            hdigit = state_diagram[sdigit, 1, cstate]

            i_bit_mask[3 * i + 2] = hdigit & 0b100
            i_bit_mask[3 * i + 1] = hdigit & 0b010
            i_bit_mask[3 * i] = hdigit & 0b001

            cstate = nstate

        # Compute ordering
        for i in range(3 * bit_length):
            order[ip] = order[ip] + i_bit_mask[i] * 2**i

    return order


def get_cpu_list(ds, X):
    """
    Return the list of the CPU intersecting with the positions
    given. Note that it will be 0-indexed.

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
    ncpu = ds.parameters["ncpu"]
    ndim = ds.parameters["ndim"]

    xmin, ymin, zmin = X.min(axis=0)
    xmax, ymax, zmax = X.max(axis=0)

    dmax = max(xmax - xmin, ymax - ymin, zmax - zmin)
    ilevel = 0
    deltax = dmax * 2

    while deltax >= dmax:
        ilevel += 1
        deltax = 0.5**ilevel

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

    idom, jdom, kdom = (np.zeros(8, dtype="int64") for _ in range(3))

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
    for i in range(ndom):
        if bit_length > 0:
            order_min = hilbert3d([idom[i], jdom[i], kdom[i]], bit_length)
        else:
            order_min = 0
        bounding_min[i] = (order_min) * dkey
        bounding_max[i] = (order_min + 1) * dkey

    bound_key = {}
    for icpu in range(1, ncpu + 1):
        bound_key[icpu - 1], bound_key[icpu] = ds.hilbert_indices[icpu]

    cpu_min, cpu_max = (np.zeros(ncpu + 1, dtype="int64") for _ in range(2))
    for icpu in range(1, ncpu + 1):
        for i in range(ndom):
            if (
                bound_key[icpu - 1] <= bounding_min[i]
                and bound_key[icpu] > bounding_min[i]
            ):
                cpu_min[i] = icpu - 1
            if (
                bound_key[icpu - 1] < bounding_max[i]
                and bound_key[icpu] >= bounding_max[i]
            ):
                cpu_max[i] = icpu

    ncpu_read = 0
    cpu_list = []
    cpu_read = np.zeros(ncpu, dtype="bool")
    for i in range(ndom):
        for j in range(cpu_min[i], cpu_max[i]):
            if not cpu_read[j]:
                ncpu_read += 1
                cpu_list.append(j)
                cpu_read[j] = True

    return sorted(cpu_list)
