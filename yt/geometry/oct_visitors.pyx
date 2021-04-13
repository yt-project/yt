# distutils: include_dirs = LIB_DIR
# distutils: libraries = STD_LIBS
"""
Oct visitor functions




"""


cimport cython
cimport numpy as np

import numpy as np

from libc.stdlib cimport free, malloc

from yt.geometry.oct_container cimport OctInfo, OctreeContainer
from yt.utilities.lib.fp_utils cimport *
from yt.utilities.lib.geometry_utils cimport encode_morton_64bit

# Now some visitor functions

cdef class OctVisitor:
    def __init__(self, OctreeContainer octree, int domain_id = -1):
        cdef int i
        self.index = 0
        self.last = -1
        self.global_index = -1
        for i in range(3):
            self.pos[i] = -1
            self.ind[i] = -1
        self.dims = 0
        self.domain = domain_id
        self.level = -1
        self.oref = octree.oref
        self.nz = (1 << (self.oref*3))

    cdef void visit(self, Oct* o, np.uint8_t selected):
        raise NotImplementedError

# This copies an integer array from the source to the destination, based on the
# selection criteria.
cdef class CopyArrayI64(OctVisitor):
    @cython.boundscheck(False)
    @cython.initializedcheck(False)
    cdef void visit(self, Oct* o, np.uint8_t selected):
        # We should always have global_index less than our source.
        # "last" here tells us the dimensionality of the array.
        if selected == 0: return
        # There are this many records between "octs"
        self.dest[self.index, :] = self.source[
                self.ind[2], self.ind[1], self.ind[0],
                self.global_index, :]
        self.index += 1

# This copies a floating point array from the source to the destination, based
# on the selection criteria.
cdef class CopyArrayF64(OctVisitor):
    @cython.boundscheck(False)
    @cython.initializedcheck(False)
    cdef void visit(self, Oct* o, np.uint8_t selected):
        # We should always have global_index less than our source.
        # "last" here tells us the dimensionality of the array.
        if selected == 0: return
        # There are this many records between "octs"
        self.dest[self.index, :] = self.source[
                self.ind[2], self.ind[1], self.ind[0],
                self.global_index, :]
        self.index += 1

# This copies a bit array from source to the destination, based on file_ind
cdef class CopyFileIndArrayI8(OctVisitor):
    def __init__(self, OctreeContainer octree, int domain_id = -1):
        super(CopyFileIndArrayI8, self).__init__(octree, domain_id)
        self.root = -1
    @cython.boundscheck(False)
    @cython.initializedcheck(False)
    cdef void visit(self, Oct* o, np.uint8_t selected):
        if self.level == 0:
            self.root += 1
        if self.last != o.domain_ind:
            self.last = o.domain_ind
            self.dest[o.domain_ind] = self.source[self.root]
            self.index += 1

# This counts the number of octs, selected or not, that the selector hits.
# Note that the selector will not recursively visit unselected octs, so this is
# still useful.
cdef class CountTotalOcts(OctVisitor):
    @cython.boundscheck(False)
    @cython.initializedcheck(False)
    cdef void visit(self, Oct* o, np.uint8_t selected):
        # Count even if not selected.
        # Number of *octs* visited.
        if self.last != o.domain_ind:
            self.index += 1
            self.last = o.domain_ind

# This counts the number of selected cells.
cdef class CountTotalCells(OctVisitor):
    @cython.boundscheck(False)
    @cython.initializedcheck(False)
    cdef void visit(self, Oct* o, np.uint8_t selected):
        # Number of *cells* visited and selected.
        self.index += selected

# Every time a cell is visited, mark it.  This will be for all visited octs.
cdef class MarkOcts(OctVisitor):
    @cython.boundscheck(False)
    @cython.initializedcheck(False)
    cdef void visit(self, Oct* o, np.uint8_t selected):
        # We mark them even if they are not selected
        if self.last != o.domain_ind:
            self.last = o.domain_ind
            self.index += 1
        self.mark[self.index, self.ind[2], self.ind[1], self.ind[0]] = 1

# Mask all the selected cells.
cdef class MaskOcts(OctVisitor):
    @cython.boundscheck(False)
    @cython.initializedcheck(False)
    cdef void visit(self, Oct* o, np.uint8_t selected):
        if selected == 0: return
        self.mask[self.global_index, self.ind[2], self.ind[1], self.ind[0]] = 1

# Compute a mapping from domain_ind to flattened index.
cdef class IndexOcts(OctVisitor):
    @cython.boundscheck(False)
    @cython.initializedcheck(False)
    cdef void visit(self, Oct* o, np.uint8_t selected):
        # Note that we provide an index even if the cell is not selected.
        if self.last != o.domain_ind:
            self.last = o.domain_ind
            self.oct_index[o.domain_ind] = self.index
            self.index += 1

# Compute a mapping from domain_ind to flattend index with some octs masked.
cdef class MaskedIndexOcts(OctVisitor):
    @cython.boundscheck(False)
    @cython.initializedcheck(False)
    cdef void visit(self, Oct* o, np.uint8_t selected):
        # Note that we provide an index even if the cell is not selected.
        if self.last != o.domain_ind:
            self.last = o.domain_ind
            if self.oct_mask[o.domain_ind] == 1:
                self.oct_index[o.domain_ind] = self.index
                self.index += 1

# Compute a mapping from domain_ind to flattened index checking mask.
cdef class IndexMaskMapOcts(OctVisitor):
    def __init__(self, OctreeContainer octree, int domain_id = -1):
        super(IndexMaskMapOcts, self).__init__(octree, domain_id)
        self.map_index = 0
    @cython.boundscheck(False)
    @cython.initializedcheck(False)
    cdef void visit(self, Oct* o, np.uint8_t selected):
        if self.last != o.domain_ind:
            self.last = o.domain_ind
            if self.oct_mask[o.domain_ind] == 1:
                if self.map_domain_ind[self.map_index] >= 0:
                    self.oct_index[self.map_domain_ind[self.map_index]] = self.index
                self.map_index += 1
            self.index += 1

# Integer coordinates
cdef class ICoordsOcts(OctVisitor):
    @cython.boundscheck(False)
    @cython.initializedcheck(False)
    cdef void visit(self, Oct* o, np.uint8_t selected):
        if selected == 0: return
        cdef int i
        for i in range(3):
            self.icoords[self.index,i] = (self.pos[i] << self.oref) + self.ind[i]
        self.index += 1

# Level
cdef class IResOcts(OctVisitor):
    @cython.boundscheck(False)
    @cython.initializedcheck(False)
    cdef void visit(self, Oct* o, np.uint8_t selected):
        if selected == 0: return
        self.ires[self.index] = self.level
        self.index += 1

# Floating point coordinates
cdef class FCoordsOcts(OctVisitor):
    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.initializedcheck(False)
    cdef void visit(self, Oct* o, np.uint8_t selected):
        # Note that this does not actually give the correct floating point
        # coordinates.  It gives them in some unit system where the domain is 1.0
        # in all directions, and assumes that they will be scaled later.
        if selected == 0: return
        cdef int i
        cdef np.float64_t c, dx
        dx = 1.0 / ((1 << self.oref) << self.level)
        for i in range(3):
            c = <np.float64_t> ((self.pos[i] << self.oref ) + self.ind[i])
            self.fcoords[self.index,i] = (c + 0.5) * dx
        self.index += 1

# Floating point widths; domain modifications are done later.
cdef class FWidthOcts(OctVisitor):
    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.initializedcheck(False)
    cdef void visit(self, Oct* o, np.uint8_t selected):
        # Note that this does not actually give the correct floating point
        # coordinates.  It gives them in some unit system where the domain is 1.0
        # in all directions, and assumes that they will be scaled later.
        if selected == 0: return
        cdef int i
        cdef np.float64_t dx
        dx = 1.0 / ((1 << self.oref) << self.level)
        for i in range(3):
            self.fwidth[self.index,i] = dx
        self.index += 1

# Mark which domains are touched by a selector.
cdef class IdentifyOcts(OctVisitor):
    @cython.boundscheck(False)
    @cython.initializedcheck(False)
    cdef void visit(self, Oct* o, np.uint8_t selected):
        # We assume that our domain has *already* been selected by, which means
        # we'll get all cells within the domain for a by-domain selector and all
        # cells within the domain *and* selector for the selector itself.
        if selected == 0: return
        self.domain_mask[o.domain - 1] = 1

# Assign domain indices to octs
cdef class AssignDomainInd(OctVisitor):
    @cython.boundscheck(False)
    @cython.initializedcheck(False)
    cdef void visit(self, Oct* o, np.uint8_t selected):
        o.domain_ind = self.global_index
        self.index += 1

# From the file, fill in C order
cdef class FillFileIndicesO(OctVisitor):
    @cython.boundscheck(False)
    @cython.initializedcheck(False)
    cdef void visit(self, Oct* o, np.uint8_t selected):
        # We fill these arrays, then inside the level filler we use these as
        # indices as we fill a second array from the self.
        if selected == 0: return
        self.levels[self.index] = self.level
        self.file_inds[self.index] = o.file_ind
        self.cell_inds[self.index] = self.oind()
        self.index +=1

# From the file, fill in F order
cdef class FillFileIndicesR(OctVisitor):
    @cython.boundscheck(False)
    @cython.initializedcheck(False)
    cdef void visit(self, Oct* o, np.uint8_t selected):
        # We fill these arrays, then inside the level filler we use these as
        # indices as we fill a second array from the self.
        if selected == 0: return
        self.levels[self.index] = self.level
        self.file_inds[self.index] = o.file_ind
        self.cell_inds[self.index] = self.rind()
        self.index +=1

# Count octs by domain
cdef class CountByDomain(OctVisitor):
    @cython.boundscheck(False)
    @cython.initializedcheck(False)
    cdef void visit(self, Oct* o, np.uint8_t selected):
        if selected == 0: return
        # NOTE: We do this for every *cell*.
        self.domain_counts[o.domain - 1] += 1

# Store the refinement mapping of the octree to be loaded later
cdef class StoreOctree(OctVisitor):
    @cython.boundscheck(False)
    @cython.initializedcheck(False)
    cdef void visit(self, Oct* o, np.uint8_t selected):
        cdef np.uint8_t res, ii
        ii = cind(self.ind[0], self.ind[1], self.ind[2])
        if o.children == NULL:
            # Not refined.
            res = 0
        else:
            res = 1
        self.ref_mask[self.index] = res
        self.index += 1

# Go from a refinement mapping to a new octree
cdef class LoadOctree(OctVisitor):
    @cython.boundscheck(False)
    @cython.initializedcheck(False)
    cdef void visit(self, Oct* o, np.uint8_t selected):
        cdef int i, ii
        ii = cind(self.ind[0], self.ind[1], self.ind[2])
        if self.ref_mask[self.index] == 0:
            # We only want to do this once.  Otherwise we end up with way too many
            # nfinest for our tastes.
            if o.file_ind == -1:
                o.children = NULL
                o.file_ind = self.nfinest[0]
                o.domain = 1
                self.nfinest[0] += 1
        elif self.ref_mask[self.index] > 0:
            if self.ref_mask[self.index] != 1 and self.ref_mask[self.index] != 8:
                print("ARRAY CLUE: ", self.ref_mask[self.index], "UNKNOWN")
                raise RuntimeError
            if o.children == NULL:
                o.children = <Oct **> malloc(sizeof(Oct *) * 8)
                for i in range(8):
                    o.children[i] = NULL
            for i in range(8):
                o.children[ii + i] = &self.octs[self.nocts[0]]
                o.children[ii + i].domain_ind = self.nocts[0]
                o.children[ii + i].file_ind = -1
                o.children[ii + i].domain = -1
                o.children[ii + i].children = NULL
                self.nocts[0] += 1
        else:
            print("SOMETHING IS AMISS", self.index)
            raise RuntimeError
        self.index += 1

cdef class MortonIndexOcts(OctVisitor):
    @cython.boundscheck(False)
    @cython.initializedcheck(False)
    cdef void visit(self, Oct* o, np.uint8_t selected):
        if selected == 0: return
        cdef np.int64_t coord[3]
        cdef int i
        for i in range(3):
            coord[i] = (self.pos[i] << self.oref) + self.ind[i]
            if (coord[i] < 0):
                raise RuntimeError("Oct coordinate in dimension {} is ".format(i)+
                                   "negative. ({})".format(coord[i]))
        self.level_arr[self.index] = self.level
        self.morton_ind[self.index] = encode_morton_64bit(
                np.uint64(coord[0]),
                np.uint64(coord[1]),
                np.uint64(coord[2]))
        self.index += 1

# Store cell index
cdef class StoreIndex(OctVisitor):
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.initializedcheck(False)
    cdef void visit(self, Oct* o, np.uint8_t selected):
        if not selected: return
        self.cell_inds[o.domain_ind, self.ind[2], self.ind[1], self.ind[0]] = self.index

        self.index += 1

cdef class BaseNeighbourVisitor(OctVisitor):
    def __init__(self, OctreeContainer octree, int domain_id=-1, int n_ghost_zones=1):
        self.octree = octree
        self.neigh_ind[0] = 0
        self.neigh_ind[1] = 0
        self.neigh_ind[2] = 0
        self.n_ghost_zones = n_ghost_zones
        super(BaseNeighbourVisitor, self).__init__(octree, domain_id)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.initializedcheck(False)
    cdef void set_neighbour_info(self, Oct *o, int ishift[3]):
        cdef int i
        cdef np.float64_t c, dx
        cdef np.int64_t ipos
        cdef np.float64_t fcoords[3]
        cdef Oct *neighbour
        cdef bint local_oct
        cdef bint other_oct
        dx = 1.0 / ((1 << self.oref) << self.level)
        local_oct = True

        # Compute position of neighbouring cell
        for i in range(3):
            c = <np.float64_t> (self.pos[i] << self.oref)
            fcoords[i] = (c + 0.5 + ishift[i]) * dx / self.octree.nn[i]
            # Assuming periodicity
            if fcoords[i] < 0:
                fcoords[i] += 1
            elif fcoords[i] > 1:
                fcoords[i] -= 1
            local_oct &= (0 <= ishift[i] <= 1)
        other_oct = not local_oct

        # Use octree to find neighbour
        if other_oct:
            neighbour = self.octree.get(fcoords, &self.oi, max_level=self.level)
        else:
            neighbour = o
            self.oi.level = self.level
            for i in range(3):
                self.oi.ipos[i] = (self.pos[i] << self.oref) + ishift[i]

        # Extra step - compute cell position in neighbouring oct (and store in oi.ipos)
        if self.oi.level == self.level - 1:
            for i in range(3):
                ipos = (((self.pos[i] << self.oref) + ishift[i])) >> 1
                if (self.oi.ipos[i] << 1) == ipos:
                    self.oi.ipos[i] = 0
                else:
                    self.oi.ipos[i] = 1
        self.neighbour = neighbour

        # Index of neighbouring cell within its oct
        for i in range(3):
            self.neigh_ind[i] = <np.uint8_t>(ishift[i]) % 2

        self.other_oct = other_oct
        if other_oct:
            if self.neighbour != NULL:
                if self.oi.level == self.level - 1:
                    # Position within neighbouring oct is stored in oi.ipos
                    for i in range(3):
                        self.neigh_ind[i] = self.oi.ipos[i]
                elif self.oi.level != self.level:
                    print('This should not happen! %s %s' % (self.oi.level, self.level))
                    self.neighbour = NULL

# Store neighbouring cell index in current cell
cdef class NeighbourCellIndexVisitor(BaseNeighbourVisitor):
    # This piece of code is very much optimizable. Here are possible routes to achieve
    # much better performance:

    # - Work oct by oct, which would reduce the number of neighbor lookup
    #   from 4³=64 to 3³=27,
    # - Use faster neighbor lookup method(s). For now, all searches are started from
    #   the root mesh down to leaf nodes, but we could instead go up the tree from the
    #   central oct then down to find all neighbors (see e.g.
    #   https://geidav.wordpress.com/2017/12/02/advanced-octrees-4-finding-neighbor-nodes/).
    # -  Pre-compute the face-neighbors of all octs.

    # Note that for the last point, algorithms exist that generate the neighbors of all
    # octs in O(1) time (https://link.springer.com/article/10.1007/s13319-015-0060-9)
    # during the octree construction.
    # Another possible solution would be to keep a unordered hash map of all the octs
    # indexed by their (3-integers) position. With such structure, finding a neighbor
    # takes O(1) time. This could even come as a replacement of the current
    # pointer-based octree structure.

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.initializedcheck(False)
    cdef void visit(self, Oct* o, np.uint8_t selected):
        cdef int i, j, k
        cdef int ishift[3]
        cdef np.uint8_t neigh_cell_ind
        cdef np.int64_t neigh_domain_ind
        if selected == 0: return
        # Work at oct level
        if self.last == o.domain_ind: return

        self.last = o.domain_ind

        cdef int i0, i1
        i0 = -self.n_ghost_zones
        i1 = 2 + self.n_ghost_zones
        # Loop over cells in and directly around oct
        for i in range(i0, i1):
            ishift[0] = i
            for j in range(i0, i1):
                ishift[1] = j
                for k in range(i0, i1):
                    ishift[2] = k
                    self.set_neighbour_info(o, ishift)

                    if not self.other_oct:
                        neigh_domain_ind = o.domain_ind
                        neigh_cell_ind   = self.neighbour_rind()
                    elif self.neighbour != NULL:
                        neigh_domain_ind = self.neighbour.domain_ind
                        neigh_cell_ind   = self.neighbour_rind()
                    else:
                        neigh_domain_ind = -1
                        neigh_cell_ind   = 8

                    self.cell_inds[self.index]   = neigh_cell_ind
                    self.domain_inds[self.index] = neigh_domain_ind

                    self.index += 1

# Store file position + cell of neighbouring cell in current cell
cdef class NeighbourCellVisitor(BaseNeighbourVisitor):

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.initializedcheck(False)
    cdef void visit(self, Oct* o, np.uint8_t selected):
        cdef int i, j, k
        cdef int ishift[3]
        cdef np.int64_t neigh_file_ind
        cdef np.uint8_t neigh_cell_ind
        cdef np.int32_t neigh_domain
        cdef np.uint8_t neigh_level
        if selected == 0: return
        # Work at oct level
        if self.last == o.domain_ind: return

        self.last = o.domain_ind

        cdef int i0, i1
        i0 = -self.n_ghost_zones
        i1 = 2 + self.n_ghost_zones
        # Loop over cells in and directly around oct
        for i in range(i0, i1):
            ishift[0] = i
            for j in range(i0, i1):
                ishift[1] = j
                for k in range(i0, i1):
                    ishift[2] = k
                    self.set_neighbour_info(o, ishift)

                    if not self.other_oct:
                        neigh_level    = self.level
                        neigh_domain   = o.domain
                        neigh_file_ind = o.file_ind
                        neigh_cell_ind = self.neighbour_rind()
                    elif self.neighbour != NULL:
                        neigh_level    = self.oi.level
                        neigh_domain   = self.neighbour.domain
                        neigh_file_ind = self.neighbour.file_ind
                        neigh_cell_ind = self.neighbour_rind()
                    else:
                        neigh_level    = 255
                        neigh_domain   = -1
                        neigh_file_ind = -1
                        neigh_cell_ind = 8

                    self.levels[self.index]    = neigh_level
                    self.file_inds[self.index] = neigh_file_ind
                    self.cell_inds[self.index] = neigh_cell_ind
                    self.domains[self.index]   = neigh_domain

                    self.index += 1
