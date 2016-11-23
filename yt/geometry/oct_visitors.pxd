"""
Oct visitor definitions file




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

cimport numpy as np

cdef struct Oct
cdef struct Oct:
    np.int64_t file_ind     # index with respect to the order in which it was
                            # added
    np.int64_t domain_ind   # index within the global set of domains
    np.int64_t domain       # (opt) addl int index
    Oct **children          # Up to 8 long

cdef struct OctPadded:
    np.int64_t file_ind
    np.int64_t domain_ind
    np.int64_t domain
    np.int64_t padding

cdef class OctVisitor:
    cdef np.uint64_t index
    cdef np.uint64_t last
    cdef np.int64_t global_index
    cdef np.int64_t pos[3]       # position in ints
    cdef np.uint8_t ind[3]              # cell position
    cdef int dims
    cdef np.int32_t domain
    cdef np.int8_t level
    cdef np.int8_t oref # This is the level of overref.  1 => 8 zones, 2 => 64, etc.
                        # To calculate nzones, 1 << (oref * 3)
    cdef np.int32_t nz

    # There will also be overrides for the memoryviews associated with the
    # specific instance.

    cdef void visit(self, Oct*, np.uint8_t selected)

    cdef inline int oind(self):
        cdef int d = (1 << self.oref)
        return (((self.ind[0]*d)+self.ind[1])*d+self.ind[2])

    cdef inline int rind(self):
        cdef int d = (1 << self.oref)
        return (((self.ind[2]*d)+self.ind[1])*d+self.ind[0])

cdef class CountTotalOcts(OctVisitor):
    pass

cdef class CountTotalCells(OctVisitor):
    pass

cdef class MarkOcts(OctVisitor):
    # Unused
    cdef np.uint8_t[:,:,:,:] mark

cdef class MaskOcts(OctVisitor):
    cdef np.uint8_t[:,:,:,:] mask

cdef class IndexOcts(OctVisitor):
    cdef np.int64_t[:] oct_index

cdef class ICoordsOcts(OctVisitor):
    cdef np.int64_t[:,:] icoords

cdef class IResOcts(OctVisitor):
    cdef np.int64_t[:] ires

cdef class FCoordsOcts(OctVisitor):
    cdef np.float64_t[:,:] fcoords

cdef class FWidthOcts(OctVisitor):
    cdef np.float64_t[:,:] fwidth

cdef class CopyArrayI64(OctVisitor):
    cdef np.int64_t[:,:,:,:,:,:] source
    cdef np.int64_t[:,:] dest

cdef class CopyArrayF64(OctVisitor):
    cdef np.float64_t[:,:,:,:,:] source
    cdef np.float64_t[:,:] dest

cdef class IdentifyOcts(OctVisitor):
    cdef np.uint8_t[:] domain_mask

cdef class AssignDomainInd(OctVisitor):
    pass

cdef class FillFileIndicesO(OctVisitor):
    cdef np.uint8_t[:] levels
    cdef np.int64_t[:] file_inds
    cdef np.uint8_t[:] cell_inds

cdef class FillFileIndicesR(OctVisitor):
    cdef np.uint8_t[:] levels
    cdef np.int64_t[:] file_inds
    cdef np.uint8_t[:] cell_inds

cdef class CountByDomain(OctVisitor):
    cdef np.int64_t[:] domain_counts

cdef class StoreOctree(OctVisitor):
    cdef np.uint8_t[:] ref_mask

cdef class LoadOctree(OctVisitor):
    cdef np.uint8_t[:] ref_mask
    cdef Oct* octs
    cdef np.int64_t *nocts
    cdef np.int64_t *nfinest

cdef inline int cind(int i, int j, int k):
    # THIS ONLY WORKS FOR CHILDREN.  It is not general for zones.
    return (((i*2)+j)*2+k)

