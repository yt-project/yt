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

cdef class OctreeContainer

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
    np.uint64_t index
    np.uint64_t last
    np.int64_t global_index
    np.int64_t pos[3]       # position in ints
    np.uint8_t ind[3]              # cell position
    int dims
    np.int32_t domain
    np.int8_t level
    np.int8_t oref # This is the level of overref.  1 => 8 zones, 2 => 64, etc.
                   # To calculate nzones, 1 << (oref * 3)
    np.int32_t nz

    # There will also be overrides for the memoryviews associated with the
    # specific instance.

    cdef __init__(self, OctreeContainer octree)
                            
    cdef void visit(self, Oct*, np.uint8_t selected)

    cdef inline int oind(self):
        cdef int d = (1 << self.oref)
        return (((self.ind[0]*d)+self.ind[1])*d+self.ind[2])

    cdef inline int rind(self):
        cdef int d = (1 << self.oref)
        return (((self.ind[2]*d)+self.ind[1])*d+self.ind[0])

cdef class CountTotalOcts(OctVisitor)

cdef class CountTotalCells(OctVisitor)

cdef class MarkOcts(OctVisitor):
    # Unused
    np.uint8_t[:,:,:,:] mark

cdef class MaskOcts(OctVisitor):
    np.uint8_t[:,:,:,:] mask

cdef class IndexOcts(OctVisitor):
    np.int64_t[:] oct_index

cdef class ICoordsOcts(OctVisitor):
    np.int64_t[:,3] icoords

cdef class IResOcts(OctVisitor):
    np.int64_t[:,3] ires

cdef class FCoordsOcts(OctVisitor):
    np.float64_t[:,3] fcoords

cdef class FWidthOcts(OctVisitor):
    np.float64_t[:,3] fwidth

cdef fused numpy_dt:
    np.float32_t
    np.float64_t
    np.int32_t
    np.int64_t

cdef class CopyArray[numpy_dt](OctVisitor):
    numpy_dt[:,:] source
    numpy_dt[:,:] dest

cdef class IdentifyOcts(OctVisitor):
    np.uint64_t[:] domain_mask

cdef class AssignDomainInd(OctVisitor):
    pass

cdef class FillFileIndicesO(OctVisitor):
    np.uint8_t[:] levels
    np.uint8_t[:] file_inds
    np.uint8_t[:] cell_inds

cdef class FillFileIndicesR(OctVisitor):
    np.uint8_t[:] levels
    np.int64_t[:] file_inds
    np.uint8_t[:] cell_inds

cdef class CountByDomain(OctVisitor):
    np.int64_t[:] domain_counts

cdef class StoreOctree(OctVisitor):
    np.uint8_t[:] ref_mask

cdef class LoadOctree(OctVisitor):
    np.uint8_t[:] ref_mask
    Oct[:] octs
    np.int64_t nocts
    np.int64_t nfinest

cdef inline int cind(int i, int j, int k):
    # THIS ONLY WORKS FOR CHILDREN.  It is not general for zones.
    return (((i*2)+j)*2+k)

