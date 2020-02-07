import numpy as np
from .io_utils import hilbert3d

def get_cpu_list(ds, X):
    '''
    Return the list of the CPU intersecting with the positions
    given. Note that it will be 0-indexed.

    Parameters
    ----------
    * ds: Dataset
      The dataset containing the information
    * X: (N, ndim) float array
      An array containing positions. They should be between 0 and 1.
    '''
    X = np.atleast_2d(X)
    if X.shape[1] != 3:
        raise NotImplementedError('This function is only implemented in 3D.')

    levelmax = ds.parameters['levelmax']
    ncpu = ds.parameters['ncpu']
    ndim = ds.parameters['ndim']

    xmin, ymin, zmin = X.min(axis=0)
    xmax, ymax, zmax = X.max(axis=0)

    dmax = max(xmax-xmin, ymax-ymin, zmax-zmin)
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


    dkey = (2**(levelmax+1) / maxdom)**ndim
    ndom = 1
    if (bit_length > 0): ndom = 8

    idom, jdom, kdom = [np.zeros(8, dtype=int) for _ in range(3)]

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
            order_min = hilbert3d(np.asarray([[idom[i], jdom[i], kdom[i]]]), bit_length)
        else:
            order_min = 0
        bounding_min[i] = (order_min  )*dkey
        bounding_max[i] = (order_min+1)*dkey

    bound_key = {}
    for icpu in range(1, ncpu+1):
        bound_key[icpu-1], bound_key[icpu] = ds.hilbert_indices[icpu]

    cpu_min, cpu_max = [np.zeros(ncpu + 1, dtype=np.int) for _ in range(2)]
    for icpu in range(1, ncpu+1):
        for i in range(ndom):
            if (bound_key[icpu - 1] <= bounding_min[i] and
                bound_key[icpu    ] >  bounding_min[i]):
                cpu_min[i] = icpu-1
            if (bound_key[icpu - 1] <  bounding_max[i] and
                bound_key[icpu    ] >= bounding_max[i]):
                cpu_max[i] = icpu

    ncpu_read = 0
    cpu_list = []
    cpu_read = np.zeros(ncpu, dtype=np.bool)
    for i in range(ndom):
        for j in range(cpu_min[i], cpu_max[i]):
            if not cpu_read[j]:
                ncpu_read += 1
                cpu_list.append(j)
                cpu_read[j] = True

    return sorted(cpu_list)
