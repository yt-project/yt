"""
Particle Deposition onto Cells




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

cimport numpy as np
import numpy as np
from libc.stdlib cimport malloc, free
cimport cython
from libc.math cimport sqrt
from cpython cimport PyObject
from yt.utilities.lib.fp_utils cimport *

from oct_container cimport Oct, OctAllocationContainer, \
    OctreeContainer, OctInfo
from cpython.array cimport array, clone
from cython.view cimport memoryview as cymemview
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
        self.nvals = nvals
        self.update_values = 0 # This is the default
        self.sph_kernel = get_kernel_func(kernel_name)

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
                     int domain_offset = 0):
        cdef int nf, i, j
        if fields is None:
            fields = []
        nf = len(fields)
        cdef np.float64_t[::cython.view.indirect, ::1] field_pointers 
        if nf > 0: field_pointers = OnceIndirect(fields)
        cdef np.float64_t pos[3]
        cdef np.float64_t[:] field_vals = np.empty(nf, dtype="float64")
        cdef int dims[3]
        dims[0] = dims[1] = dims[2] = (1 << octree.oref)
        cdef int nz = dims[0] * dims[1] * dims[2]
        cdef OctInfo oi
        cdef np.int64_t offset, moff
        cdef Oct *oct
        cdef np.int64_t numpart = positions.shape[0]
        moff = octree.get_domain_offset(domain_id + domain_offset)
        for i in range(positions.shape[0]):
            # We should check if particle remains inside the Oct here
            for j in range(nf):
                field_vals[j] = field_pointers[j,i]
            for j in range(3):
                pos[j] = positions[i, j]
            # This line should be modified to have it return the index into an
            # array based on whatever cutting of the domain we have done.  This
            # may or may not include the domain indices that we have
            # previously generated.  This way we can support not knowing the
            # full octree structure.  All we *really* care about is some
            # arbitrary offset into a field value for deposition.
            oct = octree.get(pos, &oi)
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
            self.process(dims, oi.left_edge, oi.dds,
                         offset, pos, field_vals, oct.domain_ind)
            if self.update_values == 1:
                for j in range(nf):
                    field_pointers[j][i] = field_vals[j]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def process_grid(self, gobj,
                     np.ndarray[np.float64_t, ndim=2] positions,
                     fields = None):
        cdef int nf, i, j
        if fields is None:
            fields = []
        nf = len(fields)
        cdef np.float64_t[:] field_vals = np.empty(nf, dtype="float64")
        cdef np.float64_t[::cython.view.indirect, ::1] field_pointers 
        if nf > 0: field_pointers = OnceIndirect(fields)
        cdef np.float64_t pos[3]
        cdef np.int64_t gid = getattr(gobj, "id", -1)
        cdef np.float64_t dds[3]
        cdef np.float64_t left_edge[3]
        cdef int dims[3]
        for i in range(3):
            dds[i] = gobj.dds[i]
            left_edge[i] = gobj.LeftEdge[i]
            dims[i] = gobj.ActiveDimensions[i]
        for i in range(positions.shape[0]):
            # Now we process
            for j in range(nf):
                field_vals[j] = field_pointers[j,i]
            for j in range(3):
                pos[j] = positions[i, j]
            self.process(dims, left_edge, dds, 0, pos, field_vals, gid)
            if self.update_values == 1:
                for j in range(nf):
                    field_pointers[j][i] = field_vals[j]

    cdef void process(self, int dim[3], np.float64_t left_edge[3],
                      np.float64_t dds[3], np.int64_t offset,
                      np.float64_t ppos[3], np.float64_t[:] fields,
                      np.int64_t domain_ind):
        raise NotImplementedError

cdef class CountParticles(ParticleDepositOperation):
    cdef np.int64_t[:,:,:,:] count
    def initialize(self):
        # Create a numpy array accessible to python
        self.count = append_axes(
            np.zeros(self.nvals, dtype="int64", order='F'), 4)

    @cython.cdivision(True)
    cdef void process(self, int dim[3],
                      np.float64_t left_edge[3],
                      np.float64_t dds[3],
                      np.int64_t offset, # offset into IO field
                      np.float64_t ppos[3], # this particle's position
                      np.float64_t[:] fields,
                      np.int64_t domain_ind
                      ):
        # here we do our thing; this is the kernel
        cdef int ii[3]
        cdef int i
        for i in range(3):
            ii[i] = <int>((ppos[i] - left_edge[i])/dds[i])
        self.count[ii[2], ii[1], ii[0], offset] += 1

    def finalize(self):
        arr = np.asarray(self.count)
        arr.shape = self.nvals
        return arr.astype("float64")

deposit_count = CountParticles

cdef class SimpleSmooth(ParticleDepositOperation):
    # Note that this does nothing at the edges.  So it will give a poor
    # estimate there, and since Octrees are mostly edges, this will be a very
    # poor SPH kernel.
    cdef np.float64_t[:,:,:,:] data
    cdef np.float64_t[:,:,:,:] temp

    def initialize(self):
        self.data = append_axes(
            np.zeros(self.nvals, dtype="float64", order='F'), 4)
        self.temp = append_axes(
            np.zeros(self.nvals, dtype="float64", order='F'), 4)

    @cython.cdivision(True)
    cdef void process(self, int dim[3],
                      np.float64_t left_edge[3],
                      np.float64_t dds[3],
                      np.int64_t offset,
                      np.float64_t ppos[3],
                      np.float64_t[:] fields,
                      np.int64_t domain_ind
                      ):
        cdef int ii[3]
        cdef int ib0[3]
        cdef int ib1[3]
        cdef int i, j, k, half_len
        cdef np.float64_t idist[3]
        cdef np.float64_t kernel_sum, dist
        # Smoothing length is fields[0]
        kernel_sum = 0.0
        for i in range(3):
            ii[i] = <int>((ppos[i] - left_edge[i])/dds[i])
            half_len = <int>(fields[0]/dds[i]) + 1
            ib0[i] = ii[i] - half_len
            ib1[i] = ii[i] + half_len
            if ib0[i] >= dim[i] or ib1[i] <0:
                return
            ib0[i] = iclip(ib0[i], 0, dim[i] - 1)
            ib1[i] = iclip(ib1[i], 0, dim[i] - 1)
        for i from ib0[0] <= i <= ib1[0]:
            idist[0] = (ii[0] - i) * dds[0]
            idist[0] *= idist[0]
            for j from ib0[1] <= j <= ib1[1]:
                idist[1] = (ii[1] - j) * dds[1]
                idist[1] *= idist[1]
                for k from ib0[2] <= k <= ib1[2]:
                    idist[2] = (ii[2] - k) * dds[2]
                    idist[2] *= idist[2]
                    dist = idist[0] + idist[1] + idist[2]
                    # Calculate distance in multiples of the smoothing length
                    dist = sqrt(dist) / fields[0]
                    self.temp[k,j,i,offset] = self.sph_kernel(dist)
                    kernel_sum += self.temp[k,j,i,offset]
        # Having found the kernel, deposit accordingly into gdata
        for i from ib0[0] <= i <= ib1[0]:
            for j from ib0[1] <= j <= ib1[1]:
                for k from ib0[2] <= k <= ib1[2]:
                    dist = self.temp[k,j,i,offset] / kernel_sum
                    self.data[k,j,i,offset] += fields[1] * dist

    def finalize(self):
        return self.odata

deposit_simple_smooth = SimpleSmooth

cdef class SumParticleField(ParticleDepositOperation):
    cdef np.float64_t[:,:,:,:] sum
    def initialize(self):
        self.sum = append_axes(
            np.zeros(self.nvals, dtype="float64", order='F'), 4)

    @cython.cdivision(True)
    cdef void process(self, int dim[3],
                      np.float64_t left_edge[3],
                      np.float64_t dds[3],
                      np.int64_t offset,
                      np.float64_t ppos[3],
                      np.float64_t[:] fields,
                      np.int64_t domain_ind
                      ):
        cdef int ii[3]
        cdef int i
        for i in range(3):
            ii[i] = <int>((ppos[i] - left_edge[i]) / dds[i])
        self.sum[ii[2], ii[1], ii[0], offset] += fields[0]
        return

    def finalize(self):
        sum = np.asarray(self.sum)
        sum.shape = self.nvals
        return sum

deposit_sum = SumParticleField

cdef class StdParticleField(ParticleDepositOperation):
    # Thanks to Britton and MJ Turk for the link
    # to a single-pass STD
    # http://www.cs.berkeley.edu/~mhoemmen/cs194/Tutorials/variance.pdf
    cdef np.float64_t[:,:,:,:] mk
    cdef np.float64_t[:,:,:,:] qk
    cdef np.float64_t[:,:,:,:] i
    def initialize(self):
        # we do this in a single pass, but need two scalar
        # per cell, M_k, and Q_k and also the number of particles
        # deposited into each one
        # the M_k term
        self.mk = append_axes(
            np.zeros(self.nvals, dtype="float64", order='F'), 4)
        self.qk = append_axes(
            np.zeros(self.nvals, dtype="float64", order='F'), 4)
        self.i = append_axes(
            np.zeros(self.nvals, dtype="float64", order='F'), 4)

    @cython.cdivision(True)
    cdef void process(self, int dim[3],
                      np.float64_t left_edge[3],
                      np.float64_t dds[3],
                      np.int64_t offset,
                      np.float64_t ppos[3],
                      np.float64_t[:] fields,
                      np.int64_t domain_ind
                      ):
        cdef int ii[3]
        cdef int i, cell_index
        cdef float k, mk, qk
        for i in range(3):
            ii[i] = <int>((ppos[i] - left_edge[i])/dds[i])
        k = self.i[ii[2], ii[1], ii[0], offset]
        mk = self.mk[ii[2], ii[1], ii[0], offset]
        qk = self.qk[ii[2], ii[1], ii[0], offset]
        #print k, mk, qk, cell_index
        if k == 0.0:
            # Initialize cell values
            self.mk[ii[2], ii[1], ii[0], offset] = fields[0]
        else:
            self.mk[ii[2], ii[1], ii[0], offset] = mk + (fields[0] - mk) / k
            self.qk[ii[2], ii[1], ii[0], offset] = \
                qk + (k - 1.0) * (fields[0] - mk)**2.0 / k
        self.i[ii[2], ii[1], ii[0], offset] += 1

    def finalize(self):
        # This is the standard variance
        # if we want sample variance divide by (self.oi - 1.0)
        i = np.asarray(self.i)
        std2 = np.asarray(self.qk) / i
        std2[i == 0.0] = 0.0
        std2.shape = self.nvals
        return np.sqrt(std2)

deposit_std = StdParticleField

cdef class CICDeposit(ParticleDepositOperation):
    cdef np.float64_t[:,:,:,:] field
    cdef public object ofield
    def initialize(self):
        if not all(_ > 1 for _ in self.nvals):
            from yt.utilities.exceptions import YTBoundsDefinitionError
            raise YTBoundsDefinitionError(
                "CIC requires minimum of 2 zones in all dimensions",
                self.nvals)
        self.field = append_axes(
            np.zeros(self.nvals, dtype="float64", order='F'), 4)

    @cython.cdivision(True)
    cdef void process(self, int dim[3],
                      np.float64_t left_edge[3],
                      np.float64_t dds[3],
                      np.int64_t offset, # offset into IO field
                      np.float64_t ppos[3], # this particle's position
                      np.float64_t[:] fields,
                      np.int64_t domain_ind
                      ):

        cdef int i, j, k
        cdef np.uint64_t ii
        cdef int ind[3]
        cdef np.float64_t rpos[3]
        cdef np.float64_t rdds[3][2]
        cdef np.float64_t fact, edge0, edge1, edge2
        cdef np.float64_t le0, le1, le2
        cdef np.float64_t dx, dy, dz, dx2, dy2, dz2

        # Compute the position of the central cell
        for i in range(3):
            rpos[i] = (ppos[i]-left_edge[i])/dds[i]
            rpos[i] = fclip(rpos[i], 0.5001, dim[i]-0.5001)
            ind[i] = <int> (rpos[i] + 0.5)
            # Note these are 1, then 0
            rdds[i][1] = (<np.float64_t> ind[i]) + 0.5 - rpos[i]
            rdds[i][0] = 1.0 - rdds[i][1]

        for i in range(2):
            for j in range(2):
                for k in range(2):
                    self.field[ind[2] - k, ind[1] - j, ind[0] - i, offset] += \
                        fields[0]*rdds[0][i]*rdds[1][j]*rdds[2][k]

    def finalize(self):
        rv = np.asarray(self.field)
        rv.shape = self.nvals
        return rv

deposit_cic = CICDeposit

cdef class WeightedMeanParticleField(ParticleDepositOperation):
    # Deposit both mass * field and mass into two scalars
    # then in finalize divide mass * field / mass
    cdef np.float64_t[:,:,:,:] wf
    cdef np.float64_t[:,:,:,:] w
    def initialize(self):
        self.wf = append_axes(
            np.zeros(self.nvals, dtype='float64', order='F'), 4)
        self.w = append_axes(
            np.zeros(self.nvals, dtype='float64', order='F'), 4)

    @cython.cdivision(True)
    cdef void process(self, int dim[3],
                      np.float64_t left_edge[3],
                      np.float64_t dds[3],
                      np.int64_t offset,
                      np.float64_t ppos[3],
                      np.float64_t[:] fields,
                      np.int64_t domain_ind
                      ):
        cdef int ii[3]
        cdef int i
        for i in range(3):
            ii[i] = <int>((ppos[i] - left_edge[i]) / dds[i])
        self.w[ii[2], ii[1], ii[0], offset] += fields[1]
        self.wf[ii[2], ii[1], ii[0], offset] += fields[0] * fields[1]

    def finalize(self):
        wf = np.asarray(self.wf)
        w = np.asarray(self.w)
        rv = wf / w
        rv.shape = self.nvals
        return rv

deposit_weighted_mean = WeightedMeanParticleField

cdef class MeshIdentifier(ParticleDepositOperation):
    # This is a tricky one!  What it does is put into the particle array the
    # value of the oct or block (grids will always be zero) identifier that a
    # given particle resides in
    def initialize(self):
        self.update_values = 1

    @cython.cdivision(True)
    cdef void process(self, int dim[3],
                      np.float64_t left_edge[3],
                      np.float64_t dds[3],
                      np.int64_t offset,
                      np.float64_t ppos[3],
                      np.float64_t[:] fields,
                      np.int64_t domain_ind
                      ):
        fields[0] = domain_ind

    def finalize(self):
        return

deposit_mesh_id = MeshIdentifier

cdef class NNParticleField(ParticleDepositOperation):
    cdef np.float64_t[:,:,:,:] nnfield
    cdef np.float64_t[:,:,:,:] distfield
    def initialize(self):
        self.nnfield = append_axes(
            np.zeros(self.nvals, dtype="float64", order='F'), 4)
        self.distfield = append_axes(
            np.zeros(self.nvals, dtype="float64", order='F'), 4)
        self.distfield[:] = np.inf

    @cython.cdivision(True)
    cdef void process(self, int dim[3],
                      np.float64_t left_edge[3],
                      np.float64_t dds[3],
                      np.int64_t offset,
                      np.float64_t ppos[3],
                      np.float64_t[:] fields,
                      np.int64_t domain_ind
                      ):
        # This one is a bit slow.  Every grid cell is going to be iterated
        # over, and we're going to deposit particles in it.
        cdef int i, j, k
        cdef int ii[3]
        cdef np.float64_t r2
        cdef np.float64_t gpos[3]
        gpos[0] = left_edge[0] + 0.5 * dds[0]
        for i in range(dim[0]):
            gpos[1] = left_edge[1] + 0.5 * dds[1]
            for j in range(dim[1]):
                gpos[2] = left_edge[2] + 0.5 * dds[2]
                for k in range(dim[2]):
                    r2 = ((ppos[0] - gpos[0])*(ppos[0] - gpos[0]) +
                          (ppos[1] - gpos[1])*(ppos[1] - gpos[1]) +
                          (ppos[2] - gpos[2])*(ppos[2] - gpos[2]))
                    if r2 < self.distfield[k,j,i,offset]:
                        self.distfield[k,j,i,offset] = r2
                        self.nnfield[k,j,i,offset] = fields[0]
                    gpos[2] += dds[2]
                gpos[1] += dds[1]
            gpos[0] += dds[0]
        return

    def finalize(self):
        nn = np.asarray(self.nnfield)
        nn.shape = self.nvals
        return nn

deposit_nearest = NNParticleField
