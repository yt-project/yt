cimport cython
from oct_visitors cimport StoreIndex, NeighbourVisitor, NeighbourCellVisitor
from selection_routines cimport SelectorObject, AlwaysSelector, OctreeSubsetSelector
cimport numpy as np
import numpy as np


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

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def fill_level_with_domain(
                   self, int level,
                   np.uint8_t[:] levels,
                   np.uint8_t[:] cell_inds,
                   np.int64_t[:] file_inds,
                   np.int64_t[:] domains,
                   dict dest_fields, dict source_fields,
                   np.int32_t domain,
                   np.int64_t offset = 0
                   ):
        cdef np.ndarray[np.float64_t, ndim=2] source
        cdef np.ndarray[np.float64_t, ndim=1] dest
        cdef np.float64_t tmp
        cdef int i, count
        cdef str key
        for key in dest_fields:
            dest = dest_fields[key]
            source = source_fields[key]
            count = 0
            for i in range(levels.shape[0]):
                if levels[i] != level or domains[i] != domain: continue
                count += 1
                if file_inds[i] < 0:
                    dest[i + offset] = np.nan
                else:
                    # print(f'\t{i}: Accessing source {file_inds[i]}:{cell_inds[i]} source.shape=({source.shape[0]},{source.shape[1]})')
                    tmp =source[file_inds[i], cell_inds[i]]
                    dest[i + offset] = tmp # source[file_inds[i], cell_inds[i]]
        return count

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def file_index_octs_with_shift(
            self, SelectorObject selector, int domain_id,
            int num_cells = -1):
        cdef np.int64_t i
        cdef int num_octs
        if num_cells < 0:
            num_octs = selector.count_octs(self, domain_id)
            num_cells = num_octs * 4**3
        cdef NeighbourCellVisitor visitor

        cdef np.ndarray[np.uint8_t, ndim=1] shifted_levels
        cdef np.ndarray[np.uint8_t, ndim=1] shifted_cell_inds
        cdef np.ndarray[np.int64_t, ndim=1] shifted_file_inds
        cdef np.ndarray[np.int32_t, ndim=1] neigh_domain
        shifted_levels = np.full(num_cells, 255, dtype="uint8")
        shifted_file_inds = np.full(num_cells, -1, dtype="int64")
        shifted_cell_inds = np.full(num_cells, 8, dtype="uint8")
        neigh_domain = np.full(num_cells, -1, dtype="int32")

        visitor = NeighbourCellVisitor(self, -1)
        # output: level, file_ind and cell_ind of the neighbouring cells
        visitor.shifted_levels = shifted_levels
        visitor.shifted_file_inds = shifted_file_inds
        visitor.shifted_cell_inds = shifted_cell_inds
        visitor.neigh_domain = neigh_domain
        # direction to explore and extra parameters of the visitor
        visitor.octree = self
        visitor.last = -1

        # Compute indices
        self.visit_all_octs(selector, visitor)

        return shifted_levels, shifted_cell_inds, shifted_file_inds, neigh_domain