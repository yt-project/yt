cdef class SliceSelector(SelectorObject):
    cdef public int axis
    cdef public np.float64_t coord
    cdef public int ax, ay
    cdef public int reduced_dimensionality

    def __init__(self, dobj):
        self.axis = dobj.axis
        self.coord = _ensure_code(dobj.coord)
        # If we have a reduced dimensionality dataset, we want to avoid any
        # checks against it in the axes that are beyond its dimensionality.
        # This means that if we have a 2D dataset, *all* slices along z will
        # select all the zones.
        if self.axis >= dobj.ds.dimensionality:
            self.reduced_dimensionality = 1
        else:
            self.reduced_dimensionality = 0

        self.ax = (self.axis+1) % 3
        self.ay = (self.axis+2) % 3

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
        cdef np.uint64_t icoord
        cdef np.int32_t level = gobj.Level
        _ensure_code(gobj.LeftEdge)
        _ensure_code(gobj.dds)

        if level < self.min_level or level > self.max_level:
            return None
        else:
            child_mask = gobj.child_mask
            mask = np.zeros(gobj.ActiveDimensions, dtype=np.uint8)
            if level == self.max_level:
                this_level = 1
            for i in range(3):
                if i == self.axis:
                    icoord = <np.uint64_t>(
                        (self.coord - gobj.LeftEdge.d[i])/gobj.dds[i])
                    # clip coordinate to avoid seg fault below if we're
                    # exactly at a grid boundary
                    ind[i][0] = iclip(
                        icoord, 0, gobj.ActiveDimensions[i]-1)
                    ind[i][1] = ind[i][0] + 1
                else:
                    ind[i][0] = 0
                    ind[i][1] = gobj.ActiveDimensions[i]
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
        if self.reduced_dimensionality == 1:
            return 1
        if pos[self.axis] + 0.5*dds[self.axis] > self.coord \
           and pos[self.axis] - 0.5*dds[self.axis] - grid_eps <= self.coord:
            return 1
        return 0

    cdef int select_point(self, np.float64_t pos[3]) nogil:
        # two 0-volume constructs don't intersect
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_sphere(self, np.float64_t pos[3], np.float64_t radius) nogil:
        if self.reduced_dimensionality == 1:
            return 1
        cdef np.float64_t dist = self.periodic_difference(
            pos[self.axis], self.coord, self.axis)
        if dist*dist < radius*radius:
            return 1
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_bbox(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        if self.reduced_dimensionality == 1:
            return 1
        if left_edge[self.axis] - grid_eps <= self.coord < right_edge[self.axis]:
            return 1
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_bbox_edge(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        if self.reduced_dimensionality == 1:
            return 2
        if left_edge[self.axis] - grid_eps <= self.coord < right_edge[self.axis]:
            return 2 # a box with non-zero volume can't be inside a plane
        return 0

    def _hash_vals(self):
        return (("axis", self.axis),
                ("coord", self.coord))

    def _get_state_attnames(self):
        return ("axis", "coord", "ax", "ay", "reduced_dimensionality")

slice_selector = SliceSelector
