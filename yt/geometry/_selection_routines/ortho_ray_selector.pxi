cdef class OrthoRaySelector(SelectorObject):

    cdef public np.uint8_t px_ax
    cdef public np.uint8_t py_ax
    cdef public np.float64_t px
    cdef public np.float64_t py
    cdef public int axis

    def __init__(self, dobj):
        self.axis = dobj.axis
        self.px_ax = dobj.px_ax
        self.py_ax = dobj.py_ax
        self.px = dobj.px
        self.py = dobj.py

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def fill_mask(self, gobj):
        cdef np.ndarray[np.uint8_t, ndim=3] mask
        cdef np.ndarray[np.uint8_t, ndim=3, cast=True] child_mask
        cdef int i, j, k
        cdef int total = 0
        cdef int this_level = 0
        cdef int ind[3][2]
        cdef np.int32_t level = gobj.Level
        _ensure_code(gobj.LeftEdge)
        _ensure_code(gobj.RightEdge)
        _ensure_code(gobj.dds)

        if level < self.min_level or level > self.max_level:
            return None
        else:
            child_mask = gobj.child_mask
            mask = np.zeros(gobj.ActiveDimensions, dtype=np.uint8)
            if level == self.max_level:
                this_level = 1
            ind[self.axis][0] = 0
            ind[self.axis][1] = gobj.ActiveDimensions[self.axis]
            ind[self.px_ax][0] = \
                <int> ((self.px - (gobj.LeftEdge).to_ndarray()[self.px_ax]) /
                       gobj.dds[self.px_ax])
            ind[self.px_ax][1] = ind[self.px_ax][0] + 1
            ind[self.py_ax][0] = \
                <int> ((self.py - (gobj.LeftEdge).to_ndarray()[self.py_ax]) /
                       gobj.dds[self.py_ax])
            ind[self.py_ax][1] = ind[self.py_ax][0] + 1

            with nogil:
                for i in range(ind[0][0], ind[0][1]):
                    for j in range(ind[1][0], ind[1][1]):
                        for k in range(ind[2][0], ind[2][1]):
                            if this_level == 1 or child_mask[i, j, k]:
                                mask[i, j, k] = 1
                                total += 1
            if total == 0: return None
            return mask.astype("bool")

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        if self.px >= pos[self.px_ax] - 0.5*dds[self.px_ax] and \
           self.px <  pos[self.px_ax] + 0.5*dds[self.px_ax] and \
           self.py >= pos[self.py_ax] - 0.5*dds[self.py_ax] and \
           self.py <  pos[self.py_ax] + 0.5*dds[self.py_ax]:
            return 1
        return 0

    cdef int select_point(self, np.float64_t pos[3]) nogil:
        # two 0-volume constructs don't intersect
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_sphere(self, np.float64_t pos[3], np.float64_t radius) nogil:
        cdef np.float64_t dx = self.periodic_difference(
            pos[self.px_ax], self.px, self.px_ax)
        cdef np.float64_t dy = self.periodic_difference(
            pos[self.py_ax], self.py, self.py_ax)
        if dx*dx + dy*dy < radius*radius:
            return 1
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_bbox(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        if left_edge[self.px_ax] <= self.px < right_edge[self.px_ax] and \
           left_edge[self.py_ax] <= self.py < right_edge[self.py_ax] :
            return 1
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_bbox_edge(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        if left_edge[self.px_ax] <= self.px < right_edge[self.px_ax] and \
           left_edge[self.py_ax] <= self.py < right_edge[self.py_ax] :
            return 2 # a box of non-zero volume can't be inside a ray
        return 0

    def _hash_vals(self):
        return (("px_ax", self.px_ax),
                ("py_ax", self.py_ax),
                ("px", self.px),
                ("py", self.py),
                ("axis", self.axis))

    def _get_state_attnames(self):
        return ("px_ax", "py_ax", "px", "py", "axis")

ortho_ray_selector = OrthoRaySelector
