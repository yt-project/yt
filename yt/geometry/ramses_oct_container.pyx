cimport cython
from oct_visitors cimport FillFileIndicesRNeighbour, StoreIndex, NeighbourVisitor
from selection_routines cimport SelectorObject, AlwaysSelector
cimport numpy as np
import numpy as np

# cdef class FillFileIndices(oct_visitors.OctVisitor):
#     cdef np.int64_t[:,:,:,:] cell_inds
#     @cython.boundscheck(False)
#     @cython.wraparound(False)
#     @cython.initializedcheck(False)
#     cdef void visit(self, Oct* o, np.uint8_t selected):
#         if selected == 0: return
#         self.cell_inds[o.domain_ind, self.ind[2], self.ind[1], self.ind[0]] = self.index
#         self.index += 1

cdef class RAMSESOctreeContainer(SparseOctreeContainer):
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef Oct neighbour_in_direction(self, OctInfo *oi, np.int64_t *nneighbours, Oct *o,
                                  bint periodicity[3]):
        pass

    def fill_index(self, SelectorObject selector = AlwaysSelector(None)):
        # Get the on-file index of each cell
        cdef StoreIndex visitor

        cdef np.int64_t[:, :, :, :] cell_inds,

        cell_inds = np.full((self.nocts, 2, 2, 2), -1, dtype=np.int64)

        visitor = StoreIndex(self, -1)
        visitor.cell_inds = cell_inds

        self.visit_all_octs(selector, visitor)

        return np.asarray(cell_inds)

    def neighbours_in_direction(self, int idim, int direction,
                                np.int64_t[:, :, :, :] cell_inds):
        """Return index on file of all neighbours in a given direction"""
        cdef SelectorObject always_selector = AlwaysSelector(None)

        # Store the index of the neighbour
        cdef NeighbourVisitor n_visitor
        cdef np.ndarray[np.int64_t, ndim=4] neigh_cell_inds = np.full_like(cell_inds, -1)
        n_visitor = NeighbourVisitor(self, -1)
        n_visitor.idim = idim
        n_visitor.direction = direction
        n_visitor.cell_inds = cell_inds
        n_visitor.neigh_cell_inds = neigh_cell_inds
        n_visitor.octree = self
        n_visitor.last = -1
        self.visit_all_octs(always_selector, n_visitor)

        return np.asarray(neigh_cell_inds)

    #@cython.boundscheck(False)
    @cython.wraparound(False)
    def copy_neighbour_data(self,
                            np.int64_t[:] icell, np.int64_t[:] nicell,
                            np.float64_t[:, :] input, np.float64_t[:, :] output,
                            int N,):
        cdef int i

        for i in range(N):
            if nicell[i] > -1 and icell[i] > -1:
                output[icell[i], :] = input[nicell[i], :]

    def file_index_octs_with_shift(self, SelectorObject selector, int domain_id,
                                   int idim, int direction, int num_cells = -1):
        """Return index on file of all neighbours in a given direction"""
                # We create oct arrays of the correct size
        cdef np.int64_t i
        if num_cells < 0:
            num_cells = selector.count_oct_cells(self, domain_id)

        # Fill value of each cell with its neighbouring value
        cdef FillFileIndicesRNeighbour neigh_visitor
        cdef np.ndarray[np.uint8_t, ndim=1] shifted_levels
        cdef np.ndarray[np.uint8_t, ndim=1] shifted_cell_inds
        cdef np.ndarray[np.int64_t, ndim=1] shifted_file_inds
        shifted_levels = np.zeros(num_cells, dtype="uint8")
        shifted_file_inds = np.zeros(num_cells, dtype="int64")
        shifted_cell_inds = np.zeros(num_cells, dtype="uint8")

        if self.fill_style == "r":
            neigh_visitor = FillFileIndicesRNeighbour(self, domain_id)
            # output: level, file_ind and cell_ind of the neighbouring cells
            neigh_visitor.shifted_levels = shifted_levels
            neigh_visitor.shifted_file_inds = shifted_file_inds
            neigh_visitor.shifted_cell_inds = shifted_cell_inds
            # direction to explore and extra parameters of the visitor
            neigh_visitor.idim = idim
            neigh_visitor.direction = direction
            neigh_visitor.octree = self
            neigh_visitor.last = -1
        elif self.fill_style == "o":
            raise NotImplementedError('C-style filling with spatial offset has not been implemented.')
        else:
            raise RuntimeError
        self.visit_all_octs(selector, neigh_visitor)
        return shifted_levels, shifted_cell_inds, shifted_file_inds

    def file_index_octs(self, SelectorObject selector, int domain_id,
                        num_cells = -1, spatial_offset=(0, 0, 0)):


        cdef int i, idim, direction
        cdef bint do_spatial_offset
        cdef np.ndarray[np.int64_t, ndim=4] neigh_cell_inds
        cdef np.ndarray source_shifted

        do_spatial_offset = False
        for i in range(3):
            if spatial_offset[i] == 1 or spatial_offset[i] == -1:
                idim = i
                direction = spatial_offset[i]
                if do_spatial_offset:
                    raise Exception(
                        'ERROR: You can only specify one spatial offset direction, got [%s, %s, %s]!' %
                        (spatial_offset[0], spatial_offset[1], spatial_offset[2]))
                do_spatial_offset = True

        if not do_spatial_offset:
            return super(RAMSESOctreeContainer, self).file_index_octs(
                selector, domain_id, num_cells)
        else:
            return self.file_index_octs_with_shift(
                selector, domain_id, idim, direction, num_cells)
