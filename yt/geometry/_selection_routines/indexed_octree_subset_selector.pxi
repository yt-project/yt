cdef class IndexedOctreeSubsetSelector(SelectorObject):
    # This is a numpy array, which will be a bool of ndim 1
    cdef np.uint64_t min_ind
    cdef np.uint64_t max_ind
    cdef public SelectorObject base_selector
    cdef int filter_bbox
    cdef np.float64_t DLE[3]
    cdef np.float64_t DRE[3]

    def __init__(self, dobj):
        self.min_ind = dobj.min_ind
        self.max_ind = dobj.max_ind
        self.base_selector = dobj.base_selector
        self.min_level = self.base_selector.min_level
        self.max_level = self.base_selector.max_level
        self.filter_bbox = 0
        if getattr(dobj.ds, "filter_bbox", False):
            self.filter_bbox = 1
        for i in range(3):
            self.DLE[i] = dobj.ds.domain_left_edge[i]
            self.DRE[i] = dobj.ds.domain_right_edge[i]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def select_grids(self,
                     np.ndarray[np.float64_t, ndim=2] left_edges,
                     np.ndarray[np.float64_t, ndim=2] right_edges,
                     np.ndarray[np.int32_t, ndim=2] levels):
        raise RuntimeError

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_sphere(self, np.float64_t pos[3], np.float64_t radius) nogil:
        return 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        return 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_point(self, np.float64_t pos[3]) nogil:
        cdef int i
        if self.filter_bbox == 0:
            return 1
        for i in range(3):
            if pos[i] < self.DLE[i] or pos[i] > self.DRE[i]:
                return 0
        return 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_bbox(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        return self.base_selector.select_bbox(left_edge, right_edge)

    cdef int select_grid(self, np.float64_t left_edge[3],
                         np.float64_t right_edge[3], np.int32_t level,
                         Oct *o = NULL) nogil:
        # Because visitors now use select_grid, we should be explicitly
        # checking this.
        return self.base_selector.select_grid(left_edge, right_edge, level, o)

    def _hash_vals(self):
        return (hash(self.base_selector), self.min_ind, self.max_ind)

indexed_octree_subset_selector = IndexedOctreeSubsetSelector
