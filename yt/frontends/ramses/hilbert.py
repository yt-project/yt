import numpy as np
import unyt
from .io_utils import hilbert3d

def get_cpu(ds, X):
    '''
    For each position, return the CPU that contains it.

    Parameters
    ----------
    ds : Dataset
    X : (N, ndim) float array
    '''
    if isinstance(X, (unyt.array.unit_array)):
        X = X.to('unitary').value
    X = np.atleast_2d(X)

    levelmax = ds.parameters['levelmax']
    boxlen = ds.parameters['boxlen']

    # Get number of bits to encode position
    scale = boxlen  # TODO: support non periodic boundaries
    nx_loc = 1
    bscale = 2**(levelmax+1) / scale
    ncode = nx_loc*int(bscale)
    for bit_length in range(1, 32+1):
        ncode /= 2
        if ncode < 1:
            break

    if bit_length == 32:
        raise Exception('This is not supported by RAMSES.')

    iX = np.array(X*bscale, dtype=np.int64)
    order = hilbert3d(iX, bit_length)

    hilbert_keys = np.zeros(ds.parameters['ncpu']+1)
    hilbert_keys[0] = ds.hilbert_indices[1][0]
    for i in range(ds.parameters['ncpu']):
        hilbert_keys[i+1] = ds.hilbert_indices[i+1][1]

    # This is slightly costly in memory for N*ncpus, but at least it is fast
    idom = np.argmax(order[:, None] < hilbert_keys, axis=1)

    return idom


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
