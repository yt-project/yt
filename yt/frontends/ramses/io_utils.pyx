# distutils: libraries = STD_LIBS
# distutils: include_dirs = LIB_DIR
cimport cython
cimport numpy as np

import numpy as np

from yt.geometry.oct_container cimport RAMSESOctreeContainer
from yt.utilities.cython_fortran_utils cimport FortranFile

from yt.utilities.exceptions import YTIllDefinedAMRData

ctypedef np.int32_t INT32_t
ctypedef np.int64_t INT64_t
ctypedef np.float64_t DOUBLE_t

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.nonecheck(False)
def read_amr(FortranFile f, dict headers,
             np.ndarray[np.int64_t, ndim=1] ngridbound, INT64_t min_level,
             RAMSESOctreeContainer oct_handler):

    cdef INT64_t ncpu, nboundary, max_level, nlevelmax, ncpu_and_bound
    cdef DOUBLE_t nx, ny, nz
    cdef INT64_t ilevel, icpu, n, ndim, skip_len
    cdef INT32_t ng, buffer_size
    cdef np.ndarray[np.int32_t, ndim=2] numbl
    cdef np.ndarray[np.float64_t, ndim=2] pos

    ndim = headers['ndim']
    numbl = headers['numbl']
    nboundary = headers['nboundary']
    nx, ny, nz = (((i-1.0)/2.0) for i in headers['nx'])
    nlevelmax = headers['nlevelmax']
    ncpu = headers['ncpu']

    ncpu_and_bound = nboundary + ncpu

    pos = np.empty((0, 3), dtype=np.float64)
    buffer_size = 0
    # Compute number of fields to skip. This should be 31 in 3 dimensions
    skip_len = (1          # father index
                + 2*ndim   # neighbor index
                + 2**ndim  # son index
                + 2**ndim  # cpu map
                + 2**ndim  # refinement map
    )
    # Initialize values
    max_level = 0
    for ilevel in range(nlevelmax):
        for icpu in range(ncpu_and_bound):
            if icpu < ncpu:
                ng = numbl[ilevel, icpu]
            else:
                ng = ngridbound[icpu - ncpu + nboundary*ilevel]

            if ng == 0:
                continue
            # Skip grid index, 'next' and 'prev' arrays (they are used
            # to build the linked list in RAMSES)
            f.skip(3)

            # Allocate more memory if required
            if ng > buffer_size:
                pos = np.empty((ng, 3), dtype="d")
                buffer_size = ng

            pos[:ng, 0] = f.read_vector("d") - nx
            pos[:ng, 1] = f.read_vector("d") - ny
            pos[:ng, 2] = f.read_vector("d") - nz

            # Skip father, neighbor, son, cpu map and refinement map
            f.skip(skip_len)
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
@cython.nonecheck(False)
cpdef read_offset(FortranFile f, INT64_t min_level, INT64_t domain_id, INT64_t nvar, dict headers):

    cdef np.ndarray[np.int64_t, ndim=2] offset, level_count
    cdef INT64_t ndim, twotondim, nlevelmax, n_levels, nboundary, ncpu, ncpu_and_bound
    cdef INT64_t ilevel, icpu, skip_len
    cdef INT32_t file_ilevel, file_ncache

    numbl = headers['numbl']
    ndim = headers['ndim']
    nboundary = headers['nboundary']
    nlevelmax = headers['nlevelmax']
    n_levels = nlevelmax - min_level
    ncpu = headers['ncpu']

    ncpu_and_bound = nboundary + ncpu
    twotondim = 2**ndim

    skip_len = twotondim * nvar

    # It goes: level, CPU, 8-variable (1 oct)
    offset = np.full((ncpu_and_bound, n_levels), -1, dtype=np.int64)
    level_count = np.zeros((ncpu_and_bound, n_levels), dtype=np.int64)

    cdef np.int64_t[:, ::1] level_count_view = level_count
    cdef np.int64_t[:, ::1] offset_view = offset

    for ilevel in range(nlevelmax):
        for icpu in range(ncpu_and_bound):
            file_ilevel = f.read_int()
            file_ncache = f.read_int()
            if file_ncache == 0:
                continue

            if file_ilevel != ilevel+1:
                raise YTIllDefinedAMRData(
                    'Cannot read offsets in file %s. The level read '
                    'from data (%s) is not coherent with the expected (%s)',
                    f.name, file_ilevel, ilevel)

            if ilevel >= min_level:
                offset_view[icpu, ilevel - min_level] = f.tell()
                level_count_view[icpu, ilevel - min_level] = <INT64_t> file_ncache
            f.skip(skip_len)

    return offset, level_count

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.nonecheck(False)
def fill_hydro(FortranFile f,
               np.ndarray[np.int64_t, ndim=2] offsets,
               np.ndarray[np.int64_t, ndim=2] level_count,
               list cpu_enumerator,
               np.ndarray[np.uint8_t, ndim=1] levels,
               np.ndarray[np.uint8_t, ndim=1] cell_inds,
               np.ndarray[np.int64_t, ndim=1] file_inds,
               INT64_t ndim, list all_fields, list fields,
               dict tr,
               RAMSESOctreeContainer oct_handler,
               np.ndarray[np.int32_t, ndim=1] domains=np.array([], dtype='int32')):
    cdef INT64_t offset
    cdef dict tmp
    cdef str field
    cdef INT64_t twotondim
    cdef int ilevel, icpu, ifield, nfields, nlevels, nc, ncpu_selected
    cdef np.ndarray[np.uint8_t, ndim=1] mask

    twotondim = 2**ndim
    nfields = len(all_fields)
    ncpu = offsets.shape[0]
    nlevels = offsets.shape[1]
    ncpu_selected = len(cpu_enumerator)

    mask = np.array([(field in fields) for field in all_fields], dtype=np.uint8)

    # Loop over levels
    for ilevel in range(nlevels):
        # Loop over cpu domains
        for icpu in cpu_enumerator:
            nc = level_count[icpu, ilevel]
            if nc == 0:
                continue
            offset = offsets[icpu, ilevel]
            if offset == -1:
                continue
            f.seek(offset)
            tmp = {}
            # Initalize temporary data container for io
            # note: we use Fortran ordering to reflect the in-file ordering
            for field in all_fields:
                tmp[field] = np.empty((nc, twotondim), dtype="float64", order='F')

            for i in range(twotondim):
                # Read the selected fields
                for ifield in range(nfields):
                    if not mask[ifield]:
                        f.skip()
                    else:
                        tmp[all_fields[ifield]][:, i] = f.read_vector('d') # i-th cell
            if ncpu_selected > 1:
                oct_handler.fill_level_with_domain(
                    ilevel, levels, cell_inds, file_inds, domains, tr, tmp, domain=icpu+1)
            else:
                oct_handler.fill_level(
                    ilevel, levels, cell_inds, file_inds, tr, tmp)
