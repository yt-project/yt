"""
Oct visitor definitions file




"""


cimport numpy as np


cdef struct Oct
cdef struct Oct:
    np.int64_t file_ind     # index with respect to the order in which it was
                            # added
    np.int64_t domain_ind   # index within the global set of domains
    np.int64_t domain       # (opt) addl int index
    Oct **children          # Up to 8 long

cdef struct OctInfo:
    np.float64_t left_edge[3]
    np.float64_t dds[3]
    np.int64_t ipos[3]
    np.int32_t level

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
    cdef np.uint8_t ind[3]       # cell position
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

cdef class MaskedIndexOcts(OctVisitor):
    cdef np.int64_t[:] oct_index
    cdef np.uint8_t[:] oct_mask

cdef class IndexMaskMapOcts(OctVisitor):
    cdef np.int64_t[:] oct_index
    cdef np.uint8_t[:] oct_mask
    cdef np.int64_t[:] map_domain_ind
    cdef np.uint64_t map_index

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

cdef class CopyFileIndArrayI8(OctVisitor):
    cdef np.int64_t root
    cdef np.uint8_t[:] source
    cdef np.uint8_t[:] dest

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
    cdef np.uint64_t *nocts
    cdef np.uint64_t *nfinest

cdef class MortonIndexOcts(OctVisitor):
    cdef np.uint8_t[:] level_arr
    cdef np.uint64_t[:] morton_ind

cdef inline int cind(int i, int j, int k) nogil:
    # THIS ONLY WORKS FOR CHILDREN.  It is not general for zones.
    return (((i*2)+j)*2+k)

from oct_container cimport OctreeContainer


cdef class StoreIndex(OctVisitor):
    cdef np.int64_t[:,:,:,:] cell_inds

# cimport oct_container
cdef class BaseNeighbourVisitor(OctVisitor):
    cdef int idim      # 0,1,2 for x,y,z
    cdef int direction # +1 for +x, -1 for -x
    cdef np.uint8_t neigh_ind[3]
    cdef bint other_oct
    cdef Oct *neighbour
    cdef OctreeContainer octree
    cdef OctInfo oi
    cdef int n_ghost_zones

    cdef void set_neighbour_info(self, Oct *o, int ishift[3])

    cdef inline np.uint8_t neighbour_rind(self):
        cdef int d = (1 << self.oref)
        return (((self.neigh_ind[2]*d)+self.neigh_ind[1])*d+self.neigh_ind[0])

cdef class NeighbourCellIndexVisitor(BaseNeighbourVisitor):
    cdef np.uint8_t[::1] cell_inds
    cdef np.int64_t[::1] domain_inds

cdef class NeighbourCellVisitor(BaseNeighbourVisitor):
    cdef np.uint8_t[::1] levels
    cdef np.int64_t[::1] file_inds
    cdef np.uint8_t[::1] cell_inds
    cdef np.int32_t[::1] domains
