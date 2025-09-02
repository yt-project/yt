# distutils: include_dirs = LIB_DIR
# distutils: libraries = STD_LIBS
"""
Particle Deposition onto Cells




"""


cimport numpy as np

import numpy as np

cimport cython
from libc.math cimport sqrt

from yt.utilities.lib.fp_utils cimport *

from .oct_container cimport Oct, OctInfo, OctreeContainer

from yt.utilities.lib.misc_utilities import OnceIndirect


cdef append_axes(np.ndarray arr, int naxes):
    if arr.ndim == naxes:
        return arr
    # Avoid copies
    arr2 = arr.view()
    arr2.shape = arr2.shape + (1,) * (naxes - arr2.ndim)
    return arr2

cdef class ParticleDepositOperation:
    def __init__(self, nvals, kernel_name):
        # nvals is a tuple containing the active dimensions of the
        # grid to deposit onto and the number of grids,
        # (nx, ny, nz, ngrids, vector_size)
        self.nvals = nvals
        self.update_values = 0 # This is the default

    def initialize(self, *args):
        raise NotImplementedError

    def finalize(self, *args):
        raise NotImplementedError

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def process_octree(self, OctreeContainer octree,
                     np.ndarray[np.int64_t, ndim=1] dom_ind,
                     np.ndarray[np.float64_t, ndim=2] positions,
                     fields = None, int domain_id = -1,
                     int domain_offset = 0, lvlmax = None):
        cdef int nf, i, j
        if fields is None:
            fields = []
        nf = len(fields)
        nvec = len(fields[0][0])
        cdef np.float64_t[::cython.view.indirect, :, ::1] field_pointers
        if nf > 0: field_pointers = OnceIndirect(fields)
        cdef np.float64_t pos[3]
        cdef np.float64_t[:, :] field_vals = np.empty((nf, nvec), dtype="float64")
        cdef int dims[3]
        dims[0] = dims[1] = dims[2] = octree.nz
        cdef OctInfo oi
        cdef np.int64_t offset, moff
        cdef Oct *oct
        cdef np.int8_t use_lvlmax
        moff = octree.get_domain_offset(domain_id + domain_offset)
        if lvlmax is None:
            use_lvlmax = False
            lvlmax = []
        else:
            use_lvlmax = True
        cdef np.ndarray[np.int32_t, ndim=1] lvlmaxval = np.asarray(lvlmax, dtype=np.int32)

        for i in range(positions.shape[0]):
            # We should check if particle remains inside the Oct here
            for j in range(nf):
                for k in range(nvec):
                    field_vals[j][k] = field_pointers[j, i, k]
            for j in range(3):
                pos[j] = positions[i, j]
            # This line should be modified to have it return the index into an
            # array based on whatever cutting of the domain we have done.  This
            # may or may not include the domain indices that we have
            # previously generated.  This way we can support not knowing the
            # full octree structure.  All we *really* care about is some
            # arbitrary offset into a field value for deposition.
            if not use_lvlmax:
                oct = octree.get(pos, &oi)
            else:
                oct = octree.get(pos, &oi, max_level=lvlmaxval[i])
            # This next line is unfortunate.  Basically it says, sometimes we
            # might have particles that belong to octs outside our domain.
            # For the distributed-memory octrees, this will manifest as a NULL
            # oct.  For the non-distributed memory octrees, we'll simply see
            # this as a domain_id that is not the current domain id.  Note that
            # this relies on the idea that all the particles in a region are
            # all fed to sequential domain subsets, which will not be true with
            # RAMSES, where we *will* miss particles that live in ghost
            # regions on other processors.  Addressing this is on the TODO
            # list.
            if oct == NULL or (domain_id > 0 and oct.domain != domain_id):
                continue
            # Note that this has to be our local index, not our in-file index.
            offset = dom_ind[oct.domain_ind - moff]
            if offset < 0: continue
            # Check that we found the oct ...
            self.process(dims, i, oi.left_edge, oi.dds,
                         offset, pos, field_vals, oct.domain_ind)
            if self.update_values == 1:
                for j in range(nf):
                    field_pointers[j][i, :] = field_vals[j, :]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def process_grid(self, gobj,
                     np.ndarray[np.float64_t, ndim=2, cast=True] positions,
                     fields = None):
        cdef int nf, i, j
        if fields is None:
            fields = []
        if positions.shape[0] == 0: return
        nf = len(fields)
        nvec = len(fields[0])
        cdef np.float64_t[:, :] field_vals = np.empty((nf, nvec), dtype="float64")
        cdef np.float64_t[::cython.view.indirect, :, ::1] field_pointers
        if nf > 0: field_pointers = OnceIndirect(fields)
        cdef np.float64_t pos[3]
        cdef np.int64_t gid = getattr(gobj, "id", -1)
        cdef np.float64_t dds[3]
        cdef np.float64_t left_edge[3]
        cdef np.float64_t right_edge[3]
        cdef int dims[3]
        for i in range(3):
            dds[i] = gobj.dds[i]
            left_edge[i] = gobj.LeftEdge[i]
            right_edge[i] = gobj.RightEdge[i]
            dims[i] = gobj.ActiveDimensions[i]
        for i in range(positions.shape[0]):
            # Now we process
            for j in range(nf):
                field_vals[j][:] = field_pointers[j, i]
            for j in range(3):
                pos[j] = positions[i, j]
            continue_loop = False
            for j in range(3):
                if pos[j] < left_edge[j] or pos[j] > right_edge[j]:
                    continue_loop = True
            if continue_loop:
                continue
            self.process(dims, i, left_edge, dds, 0, pos, field_vals, gid)
            if self.update_values == 1:
                for j in range(nf):
                    field_pointers[j][i, :] = field_vals[j, :]

    cdef int process(self, int dim[3], int ipart, np.float64_t left_edge[3],
                     np.float64_t dds[3], np.int64_t offset,
                     np.float64_t ppos[3], np.float64_t[:, :] fields,
                     np.int64_t domain_ind) except -1 nogil:
        with gil:
            raise NotImplementedError


cdef class SumParticleField(ParticleDepositOperation):
    cdef np.float64_t[:,:,:,:,:] sum
    def initialize(self):
        self.sum = append_axes(
            np.zeros(self.nvals, dtype="float64", order='F'), 5)

    @cython.cdivision(True)
    @cython.boundscheck(False)
    cdef int process(self, int dim[3], int ipart,
                     np.float64_t left_edge[3],
                     np.float64_t dds[3],
                     np.int64_t offset,
                     np.float64_t ppos[3],
                     np.float64_t[:, :] fields,
                     np.int64_t domain_ind
                     ) except -1 nogil:
        cdef int ii[3]
        cdef int i
        for i in range(3):
            ii[i] = <int>((ppos[i] - left_edge[i]) / dds[i])
        self.sum[ii[2], ii[1], ii[0], offset, :] = fields[0, :]
        return 0

    def finalize(self):
        sum = np.asarray(self.sum)
        sum.shape = self.nvals
        return sum

deposit_sum = SumParticleField
