cimport cython
cimport oct_visitors
from selection_routines cimport SelectorObject, AlwaysSelector
cimport numpy as np
import numpy as np

cdef class FillFileIndices(oct_visitors.OctVisitor):
    cdef np.int64_t[:,:,:,:] cell_inds
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.initializedcheck(False)
    cdef void visit(self, Oct* o, np.uint8_t selected):
        if selected == 0: return
        self.cell_inds[o.domain_ind, self.ind[2], self.ind[1], self.ind[0]] = self.index
        self.index += 1
    
cdef class BaseNeighbourVisitor(oct_visitors.OctVisitor):
    cdef int idim      # 0,1,2 for x,y,z
    cdef int direction # +1 for +x, -1 for -x
    cdef np.uint8_t neigh_ind[3]
    cdef np.int64_t[:,:,:,:] cell_inds
    cdef RAMSESOctreeContainer octree
    cdef OctInfo oi
    cdef Oct *neighbour

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.initializedcheck(False)
    cdef void set_neighbour_oct(self):
        cdef int i
        cdef np.float64_t c, dx
        cdef np.int64_t ipos
        cdef np.float64_t fcoords[3]
        cdef Oct *neighbour
        dx = 1.0 / ((1 << self.oref) << self.level)
        # Compute position of neighbouring cell
        for i in range(3):
            c = <np.float64_t> ((self.pos[i] << self.oref) + self.ind[i])
            if i == self.idim:
                fcoords[i] = (c + 0.5 + self.direction) * dx / self.octree.nn[i]
            else:
                fcoords[i] = (c + 0.5) * dx / self.octree.nn[i]

        # Use octree to find neighbour
        neighbour = self.octree.get(fcoords, &self.oi, max_level=self.level)

        # Extra step - compute cell position in neighbouring oct (and store in oi.ipos)
        if self.oi.level == self.level - 1:
            for i in range(3):
                ipos = (((self.pos[i] << self.oref) + self.ind[i])) >> 1
                if i == self.idim:
                    ipos += self.direction
                if (self.oi.ipos[i] << 1) == ipos:
                    self.oi.ipos[i] = 0
                else:
                    self.oi.ipos[i] = 1
        self.neighbour = neighbour

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.initializedcheck(False)
    cdef np.int64_t get_neighbour_cell_index(self, Oct* o, np.uint8_t selected):
        cdef int i
        cdef np.int64_t cell_ind
        cdef bint other_oct  # True if the neighbouring cell lies in another oct

        # Note that we provide an index even if the cell is not selected.
        # if selected == 0: return -1
        # Index of neighbouring cell within its oct
        for i in range(3):
            if i == self.idim:
                self.neigh_ind[i] = (self.ind[i] + self.direction)
                other_oct = self.neigh_ind[i] != 0 and self.neigh_ind[i] != 1
                if other_oct:
                    self.neigh_ind[i] %= 2
            else:
                self.neigh_ind[i] = self.ind[i]

        if not other_oct:
            # Simple case: the neighbouring cell is within the oct
            cell_ind = self.cell_inds[o.domain_ind, self.neigh_ind[2], self.neigh_ind[1], self.neigh_ind[0]]
        else:
            # Complicated case: the cell is in a neighbouring oct
            if self.last != o.domain_ind:
                self.set_neighbour_oct()
                self.last = o.domain_ind

            if self.neighbour != NULL:
                if self.oi.level == self.level - 1:
                    # Position within neighbouring oct is stored in oi.ipos
                    for i in range(3):
                        self.neigh_ind[i] = self.oi.ipos[i]
                elif self.oi.level != self.level:
                    print('This should not happen! %s %s' % (self.oi.level, self.level))
                    return -1
                cell_ind = self.cell_inds[self.neighbour.domain_ind, self.neigh_ind[2], self.neigh_ind[1], self.neigh_ind[0]]
            else:
                cell_ind = -1
        return cell_ind

    cdef inline np.uint8_t neighbour_rind(self):
        cdef int d = (1 << self.oref)
        return (((self.neigh_ind[2]*d)+self.neigh_ind[1])*d+self.neigh_ind[0])

cdef class NeighbourVisitor(BaseNeighbourVisitor):
    cdef np.int64_t[:,:,:,:] neigh_cell_inds
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.initializedcheck(False)
    cdef void visit(self, Oct* o, np.uint8_t selected):
        cdef np.int64_t cell_ind
        cell_ind = self.get_neighbour_cell_index(o, selected)
        self.neigh_cell_inds[o.domain_ind, self.ind[2], self.ind[1], self.ind[0]] = cell_ind

cdef class FillFileIndicesRNeighbour(BaseNeighbourVisitor):
    cdef np.uint8_t[:] shifted_levels
    cdef np.int64_t[:] shifted_file_inds
    cdef np.uint8_t[:] shifted_cell_inds

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.initializedcheck(False)
    cdef void visit(self, Oct* o, np.uint8_t selected):
        cdef np.int64_t neigbour_cell_index
        if selected == 0: return
        # Note: only selected items have an index
        neighbour_index = self.get_neighbour_cell_index(o, selected)
        self.shifted_levels[self.index] = self.oi.level
        self.shifted_file_inds[self.index] = self.neighbour.file_ind
        self.shifted_cell_inds[self.index] = self.neighbour_rind()
        self.index += 1
        # if neighbour_index > -1:
        #     self.shifted_levels[neighbour_index] = self.level
        #     self.shifted_file_inds[neighbour_index] = o.file_ind
        #     self.shifted_cell_inds[neighbour_index] = self.rind()

cdef class RAMSESOctreeContainer(SparseOctreeContainer):
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef Oct neighbour_in_direction(self, OctInfo *oi, np.int64_t *nneighbours, Oct *o,
                                  bint periodicity[3]):
        pass

    def neighbours_in_direction(self, int idim, int direction, SelectorObject selector = AlwaysSelector(None)):
        """Return index on file of all neighbours in a given direction"""
        cdef SelectorObject always_selector = AlwaysSelector(None)

        cdef int num_cells = selector.count_oct_cells(self, -1)

        # Get the on-file index of each cell
        cdef FillFileIndices visitor
        cdef np.ndarray[np.int64_t, ndim=4] cell_inds = np.zeros((self.nocts, 2, 2, 2), dtype="int64")

        visitor = FillFileIndices(self, -1)
        visitor.cell_inds = cell_inds

        self.visit_all_octs(selector, visitor)

        # Store the index of the neighbour
        cdef NeighbourVisitor n_visitor
        cdef np.ndarray[np.int64_t, ndim=4] neigh_cell_inds = np.empty_like(cell_inds)
        n_visitor = NeighbourVisitor(self, -1)
        n_visitor.idim = idim
        n_visitor.direction = direction
        n_visitor.cell_inds = cell_inds
        n_visitor.neigh_cell_inds = neigh_cell_inds
        n_visitor.octree = self
        n_visitor.last = -1
        self.visit_all_octs(always_selector, n_visitor)

        return np.asarray(cell_inds), np.asarray(neigh_cell_inds)
