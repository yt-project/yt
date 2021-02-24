cdef class ComposeSelector(SelectorObject):
    cdef SelectorObject selector1
    cdef SelectorObject selector2

    def __init__(self, dobj, selector1, selector2):
        self.selector1 = selector1
        self.selector2 = selector2
        self.min_level = max(selector1.min_level, selector2.min_level)
        self.max_level = min(selector1.max_level, selector2.max_level)

    def select_grids(self,
                     np.ndarray[np.float64_t, ndim=2] left_edges,
                     np.ndarray[np.float64_t, ndim=2] right_edges,
                     np.ndarray[np.int32_t, ndim=2] levels):
        return np.logical_or(
                    self.selector1.select_grids(left_edges, right_edges, levels),
                    self.selector2.select_grids(left_edges, right_edges, levels))

    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        if self.selector1.select_cell(pos, dds) and \
                self.selector2.select_cell(pos, dds):
            return 1
        else:
            return 0

    cdef int select_grid(self, np.float64_t left_edge[3],
                         np.float64_t right_edge[3], np.int32_t level,
                         Oct *o = NULL) nogil:
        if self.selector1.select_grid(left_edge, right_edge, level, o) or \
                self.selector2.select_grid(left_edge, right_edge, level, o):
            return 1
        else:
            return 0

    cdef int select_point(self, np.float64_t pos[3]) nogil:
        if self.selector1.select_point(pos) and \
                self.selector2.select_point(pos):
            return 1
        else:
            return 0

    cdef int select_sphere(self, np.float64_t pos[3], np.float64_t radius) nogil:
        if self.selector1.select_sphere(pos, radius) and \
                self.selector2.select_sphere(pos, radius):
            return 1
        else:
            return 0

    cdef int select_bbox(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        if self.selector1.select_bbox(left_edge, right_edge) and \
                self.selector2.select_bbox(left_edge, right_edge):
            return 1
        else:
            return 0

    cdef int select_bbox_edge(self, np.float64_t left_edge[3],
                              np.float64_t right_edge[3]) nogil:
        cdef int rv1 = self.selector1.select_bbox_edge(left_edge, right_edge)
        if rv1 == 0: return 0
        cdef int rv2 = self.selector2.select_bbox_edge(left_edge, right_edge)
        if rv2 == 0: return 0
        return max(rv1, rv2)

    def _hash_vals(self):
        return (hash(self.selector1), hash(self.selector2))

compose_selector = ComposeSelector
