cimport cython
cimport numpy as np
import numpy as np
import yt.utilities.fortran_utils as fpu

from yt.geometry.oct_container import \
    RAMSESOctreeContainer

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def read_amr(f, headers, np.ndarray[np.int64_t, ndim=1] ngridbound,
             int min_level,
             oct_handler,
):

    cdef int ncpu, nboundary, max_level, nlevelmax, ncpu_and_bound
    cdef double nx, ny, nz
    cdef int ilevel, icpu, ng
    cdef np.ndarray[np.int32_t, ndim=2] numbl
    cdef np.ndarray[np.float64_t, ndim=2] pos

    numbl = headers['numbl']
    nboundary = headers['nboundary']
    nx, ny, nz = (((i-1.0)/2.0) for i in headers['nx'])
    nlevelmax = headers['nlevelmax']
    ncpu = headers['ncpu']

    ncpu_and_bound = nboundary + ncpu

    pos = np.empty((0, 3), dtype=np.float64)
    for ilevel in range(nlevelmax):
        for icpu in range(ncpu_and_bound):
            if icpu < ncpu:
                ng = numbl[ilevel, icpu]
            else:
                ng = ngridbound[icpu - ncpu + nboundary*ilevel]

            if ng == 0:
                continue
            # ind = fpu.read_vector(f, "I").astype("int64")  # NOQA
            fpu.skip(f, 3)

            # Reallocate memory if requred
            if ng > pos.shape[0]:
                pos = np.empty((ng, 3), dtype="d")

            pos[:ng, 0] = fpu.read_vector(f, "d") - nx
            pos[:ng, 1] = fpu.read_vector(f, "d") - ny
            pos[:ng, 2] = fpu.read_vector(f, "d") - nz

            # pos = np.empty((ng, 3), dtype='float64')
            # pos[:,0] = fpu.read_vector(f, "d") - nx
            # pos[:,1] = fpu.read_vector(f, "d") - ny
            # pos[:,2] = fpu.read_vector(f, "d") - nz

            fpu.skip(f, 31)
            # Note that we're adding *grids*, not individual cells.
            if ilevel >= min_level:
                n = oct_handler.add(icpu + 1, ilevel - min_level, pos[:ng, :],
                                    count_boundary = 1)
                if n > 0:
                    max_level = max(ilevel - min_level, max_level)

    return max_level

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def read_offset(f, int min_level, int domain_id, int nvar, headers):

    cdef np.ndarray[np.int64_t, ndim=1] offset, level_count
    cdef int ndim, twotondim, nlevelmax, n_levels, nboundary, ncpu, ncpu_and_bound
    cdef int file_ilevel, file_ncache

    numbl = headers['numbl']
    ndim = headers['ndim']
    nboundary = headers['nboundary']
    nlevelmax = headers['nlevelmax']
    n_levels = nlevelmax - min_level
    ncpu = headers['ncpu']

    ncpu_and_bound = nboundary + ncpu
    twotondim = 2**ndim

    # It goes: level, CPU, 8-variable (1 cube)
    offset = np.zeros(n_levels, dtype=np.int64)
    offset -= 1
    level_count = np.zeros(n_levels, dtype=np.int64)
    for ilevel in range(nlevelmax):
        for cpu in range(ncpu_and_bound):
            file_ilevel = fpu.read_vector(f, 'I')
            file_ncache = fpu.read_vector(f, 'I')
            if file_ncache == 0: continue
            assert(file_ilevel == ilevel+1)
            if cpu + 1 == domain_id and ilevel >= min_level:
                offset[ilevel - min_level] = f.tell()
                level_count[ilevel - min_level] = file_ncache
            fpu.skip(f, twotondim * nvar)

    return offset, level_count
