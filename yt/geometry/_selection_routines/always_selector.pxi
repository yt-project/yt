cdef class AlwaysSelector(SelectorObject):

    def __init__(self, dobj):
        self.overlap_cells = 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def select_grids(self,
                     np.ndarray[np.float64_t, ndim=2] left_edges,
                     np.ndarray[np.float64_t, ndim=2] right_edges,
                     np.ndarray[np.int32_t, ndim=2] levels):
        cdef int ng = left_edges.shape[0]
        cdef np.ndarray[np.uint8_t, ndim=1] gridi = np.ones(ng, dtype='uint8')
        return gridi.astype("bool")

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        return 1

    cdef int select_grid(self, np.float64_t left_edge[3],
                         np.float64_t right_edge[3], np.int32_t level,
                         Oct *o = NULL) nogil:
        return 1

    cdef int select_point(self, np.float64_t pos[3]) nogil:
        return 1

    cdef int select_sphere(self, np.float64_t pos[3], np.float64_t radius) nogil:
        return 1

    cdef int select_bbox(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        return 1

    cdef int select_bbox_edge(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        return 1

    def _hash_vals(self):
        return ("always", 1,)

always_selector = AlwaysSelector
