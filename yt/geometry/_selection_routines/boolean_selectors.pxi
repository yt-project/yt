
cdef class BooleanSelector(SelectorObject):

    def __init__(self, dobj):
        # Note that this has a different API than the other selector objects,
        # so will not work as a traditional data selector.
        if not hasattr(dobj.dobj1, "selector"):
            self.sel1 = dobj.dobj1
        else:
            self.sel1 = dobj.dobj1.selector
        if not hasattr(dobj.dobj2, "selector"):
            self.sel2 = dobj.dobj2
        else:
            self.sel2 = dobj.dobj2.selector

cdef class BooleanANDSelector(BooleanSelector):
    cdef int select_bbox(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        cdef int rv1 = self.sel1.select_bbox(left_edge, right_edge)
        if rv1 == 0: return 0
        cdef int rv2 = self.sel2.select_bbox(left_edge, right_edge)
        if rv2 == 0: return 0
        return 1

    cdef int select_bbox_edge(self, np.float64_t left_edge[3],
                              np.float64_t right_edge[3]) nogil:
        cdef int rv1 = self.sel1.select_bbox_edge(left_edge, right_edge)
        if rv1 == 0: return 0
        cdef int rv2 = self.sel2.select_bbox_edge(left_edge, right_edge)
        if rv2 == 0: return 0
        return max(rv1, rv2)

    cdef int select_grid(self, np.float64_t left_edge[3],
                         np.float64_t right_edge[3], np.int32_t level,
                         Oct *o = NULL) nogil:
        cdef int rv1 = self.sel1.select_grid(left_edge, right_edge, level, o)
        if rv1 == 0: return 0
        cdef int rv2 = self.sel2.select_grid(left_edge, right_edge, level, o)
        if rv2 == 0: return 0
        return 1

    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        cdef int rv1 = self.sel1.select_cell(pos, dds)
        if rv1 == 0: return 0
        cdef int rv2 = self.sel2.select_cell(pos, dds)
        if rv2 == 0: return 0
        return 1

    cdef int select_point(self, np.float64_t pos[3]) nogil:
        cdef int rv1 = self.sel1.select_point(pos)
        if rv1 == 0: return 0
        cdef int rv2 = self.sel2.select_point(pos)
        if rv2 == 0: return 0
        return 1

    cdef int select_sphere(self, np.float64_t pos[3], np.float64_t radius) nogil:
        cdef int rv1 = self.sel1.select_sphere(pos, radius)
        if rv1 == 0: return 0
        cdef int rv2 = self.sel2.select_sphere(pos, radius)
        if rv2 == 0: return 0
        return 1

    def _hash_vals(self):
        return (self.sel1._hash_vals() +
                ("and",) +
                self.sel2._hash_vals())

cdef class BooleanORSelector(BooleanSelector):
    cdef int select_bbox(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        cdef int rv1 = self.sel1.select_bbox(left_edge, right_edge)
        if rv1 == 1: return 1
        cdef int rv2 = self.sel2.select_bbox(left_edge, right_edge)
        if rv2 == 1: return 1
        return 0

    cdef int select_bbox_edge(self, np.float64_t left_edge[3],
                              np.float64_t right_edge[3]) nogil:
        cdef int rv1 = self.sel1.select_bbox_edge(left_edge, right_edge)
        if rv1 == 1: return 1
        cdef int rv2 = self.sel2.select_bbox_edge(left_edge, right_edge)
        if rv2 == 1: return 1
        return max(rv1, rv2)

    cdef int select_grid(self, np.float64_t left_edge[3],
                         np.float64_t right_edge[3], np.int32_t level,
                         Oct *o = NULL) nogil:
        cdef int rv1 = self.sel1.select_grid(left_edge, right_edge, level, o)
        if rv1 == 1: return 1
        cdef int rv2 = self.sel2.select_grid(left_edge, right_edge, level, o)
        if rv2 == 1: return 1
        if (rv1 == 2) or (rv2 == 2): return 2
        return 0

    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        cdef int rv1 = self.sel1.select_cell(pos, dds)
        if rv1 == 1: return 1
        cdef int rv2 = self.sel2.select_cell(pos, dds)
        if rv2 == 1: return 1
        return 0

    cdef int select_point(self, np.float64_t pos[3]) nogil:
        cdef int rv1 = self.sel1.select_point(pos)
        if rv1 == 1: return 1
        cdef int rv2 = self.sel2.select_point(pos)
        if rv2 == 1: return 1
        return 0

    cdef int select_sphere(self, np.float64_t pos[3], np.float64_t radius) nogil:
        cdef int rv1 = self.sel1.select_sphere(pos, radius)
        if rv1 == 1: return 1
        cdef int rv2 = self.sel2.select_sphere(pos, radius)
        if rv2 == 1: return 1
        return 0

    def _hash_vals(self):
        return (self.sel1._hash_vals() +
                ("or",) +
                self.sel2._hash_vals())

cdef class BooleanNOTSelector(BooleanSelector):
    cdef int select_bbox(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        # We always return True here, because we don't have a "fully included"
        # check anywhere else.
        return 1

    cdef int select_bbox_edge(self, np.float64_t left_edge[3],
                              np.float64_t right_edge[3]) nogil:
        cdef int rv1 = self.sel1.select_bbox_edge(left_edge, right_edge)
        if rv1 == 0: return 1
        elif rv1 == 1: return 0
        return 2

    cdef int select_grid(self, np.float64_t left_edge[3],
                         np.float64_t right_edge[3], np.int32_t level,
                         Oct *o = NULL) nogil:
        return 1

    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        cdef int rv1 = self.sel1.select_cell(pos, dds)
        if rv1 == 0: return 1
        return 0

    cdef int select_point(self, np.float64_t pos[3]) nogil:
        cdef int rv1 = self.sel1.select_point(pos)
        if rv1 == 0: return 1
        return 0

    cdef int select_sphere(self, np.float64_t pos[3], np.float64_t radius) nogil:
        cdef int rv1 = self.sel1.select_sphere(pos, radius)
        if rv1 == 0: return 1
        return 0

    def _hash_vals(self):
        return (self.sel1._hash_vals() +
                ("not",))

cdef class BooleanXORSelector(BooleanSelector):

    cdef int select_bbox(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        # We always return True here, because we don't have a "fully included"
        # check anywhere else.
        return 1

    cdef int select_bbox_edge(self, np.float64_t left_edge[3],
                              np.float64_t right_edge[3]) nogil:
        # Return 2 in cases where one or both selectors partially overlap since
        # part of the bounding box could satisfy the condition unless the
        # selectors are identical.
        cdef int rv1 = self.sel1.select_bbox_edge(left_edge, right_edge)
        cdef int rv2 = self.sel2.select_bbox_edge(left_edge, right_edge)
        if rv1 == rv2:
            if rv1 == 2:
                # If not identical, part of the bbox will be touched by one
                # selector and not the other.
                # if self.sel1 == self.sel2: return 0  # requires gil
                return 2
            return 0
        if rv1 == 0: return rv2
        if rv2 == 0: return rv1
        return 2  # part of bbox only touched by selector fully covering bbox

    cdef int select_grid(self, np.float64_t left_edge[3],
                         np.float64_t right_edge[3], np.int32_t level,
                         Oct *o = NULL) nogil:
        return 1

    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        cdef int rv1 = self.sel1.select_cell(pos, dds)
        cdef int rv2 = self.sel2.select_cell(pos, dds)
        if rv1 == rv2: return 0
        return 1

    cdef int select_point(self, np.float64_t pos[3]) nogil:
        cdef int rv1 = self.sel1.select_point(pos)
        cdef int rv2 = self.sel2.select_point(pos)
        if rv1 == rv2: return 0
        return 1

    cdef int select_sphere(self, np.float64_t pos[3], np.float64_t radius) nogil:
        cdef int rv1 = self.sel1.select_sphere(pos, radius)
        cdef int rv2 = self.sel2.select_sphere(pos, radius)
        if rv1 == rv2: return 0
        return 1

    def _hash_vals(self):
        return (self.sel1._hash_vals() +
                ("xor",) +
                self.sel2._hash_vals())

cdef class BooleanNEGSelector(BooleanSelector):

    cdef int select_bbox(self, np.float64_t left_edge[3],
                         np.float64_t right_edge[3]) nogil:
        # We always return True here, because we don't have a "fully included"
        # check anywhere else.
        return self.sel1.select_bbox(left_edge, right_edge)

    cdef int select_bbox_edge(self, np.float64_t left_edge[3],
                              np.float64_t right_edge[3]) nogil:
        cdef int rv1 = self.sel1.select_bbox_edge(left_edge, right_edge)
        if rv1 == 0: return 0
        cdef int rv2 = self.sel2.select_bbox_edge(left_edge, right_edge)
        if rv2 == 1:
            return 0
        elif rv2 == 0:
            return rv1
        # If sel2 is partial, then sel1 - sel2 will be partial as long
        # as sel1 != sel2
        # if self.sel1 == self.sel2: return 0  # requires gil
        return 2

    cdef int select_grid(self, np.float64_t left_edge[3],
                         np.float64_t right_edge[3], np.int32_t level,
                         Oct *o = NULL) nogil:
        return self.sel1.select_grid(left_edge, right_edge, level, o)

    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        cdef int rv1 = self.sel1.select_cell(pos, dds)
        if rv1 == 0: return 0
        cdef int rv2 = self.sel2.select_cell(pos, dds)
        if rv2 == 1: return 0
        return 1

    cdef int select_point(self, np.float64_t pos[3]) nogil:
        cdef int rv1 = self.sel1.select_point(pos)
        if rv1 == 0: return 0
        cdef int rv2 = self.sel2.select_point(pos)
        if rv2 == 1: return 0
        return 1

    cdef int select_sphere(self, np.float64_t pos[3], np.float64_t radius) nogil:
        cdef int rv1 = self.sel1.select_sphere(pos, radius)
        if rv1 == 0: return 0
        cdef int rv2 = self.sel2.select_sphere(pos, radius)
        if rv2 == 1: return 0
        return 1

    def _hash_vals(self):
        return (self.sel1._hash_vals() +
                ("neg",) +
                self.sel2._hash_vals())

cdef class ChainedBooleanSelector(SelectorObject):
    cdef int n_obj
    cdef np.ndarray selectors
    def __init__(self, dobj):
        # These are data objects, not selectors
        self.n_obj = len(dobj.data_objects)
        self.selectors = np.empty(self.n_obj, dtype="object")
        for i in range(self.n_obj):
            self.selectors[i] = dobj.data_objects[i].selector

cdef class ChainedBooleanANDSelector(ChainedBooleanSelector):
    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef int select_bbox(self, np.float64_t left_edge[3],
                         np.float64_t right_edge[3]) nogil:
        with gil:
            for i in range(self.n_obj):
                if (<SelectorObject>self.selectors[i]).select_bbox(
                        left_edge, right_edge) == 0:
                    return 0
        return 1

    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef int select_bbox_edge(self, np.float64_t left_edge[3],
                              np.float64_t right_edge[3]) nogil:
        cdef int selected = 1
        cdef int ret
        with gil:
            for i in range(self.n_obj):
                ret = (<SelectorObject>self.selectors[i]).select_bbox_edge(
                    left_edge, right_edge)
                if ret == 0:
                    return 0
                elif ret == 2:
                    selected = 2
        return selected

    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef int select_grid(self, np.float64_t left_edge[3],
                         np.float64_t right_edge[3], np.int32_t level,
                         Oct *o = NULL) nogil:
        with gil:
            for i in range(self.n_obj):
                if (<SelectorObject>self.selectors[i]).select_grid(
                        left_edge, right_edge, level, o) == 0:
                    return 0
        return 1

    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        with gil:
            for i in range(self.n_obj):
                if (<SelectorObject>self.selectors[i]).select_cell(
                        pos, dds) == 0:
                    return 0
        return 1

    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef int select_point(self, np.float64_t pos[3]) nogil:
        with gil:
            for i in range(self.n_obj):
                if (<SelectorObject>self.selectors[i]).select_point(pos) == 0:
                    return 0
        return 1

    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef int select_sphere(self, np.float64_t pos[3], np.float64_t radius) nogil:
        with gil:
            for i in range(self.n_obj):
                if (<SelectorObject>self.selectors[i]).select_sphere(
                        pos, radius) == 0:
                    return 0
        return 1

    def _hash_vals(self):
        v = ("chained_and",)
        for s in self.selectors:
            v += s._hash_vals()
        return v

intersection_selector = ChainedBooleanANDSelector

cdef class ChainedBooleanORSelector(ChainedBooleanSelector):
    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef int select_bbox(self, np.float64_t left_edge[3],
                         np.float64_t right_edge[3]) nogil:
        with gil:
            for i in range(self.n_obj):
                if (<SelectorObject>self.selectors[i]).select_bbox(
                        left_edge, right_edge) == 1:
                    return 1
        return 0

    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef int select_bbox_edge(self, np.float64_t left_edge[3],
                         np.float64_t right_edge[3]) nogil:
        cdef int selected = 0
        cdef int ret
        with gil:
            for i in range(self.n_obj):
                ret = (<SelectorObject>self.selectors[i]).select_bbox_edge(
                    left_edge, right_edge)
                if ret == 1:
                    return 1
                elif ret == 2:
                    selected = 2
        return selected

    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef int select_grid(self, np.float64_t left_edge[3],
                         np.float64_t right_edge[3], np.int32_t level,
                         Oct *o = NULL) nogil:
        with gil:
            for i in range(self.n_obj):
                if (<SelectorObject>self.selectors[i]).select_grid(
                        left_edge, right_edge, level, o) == 1:
                    return 1
        return 0

    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        with gil:
            for i in range(self.n_obj):
                if (<SelectorObject>self.selectors[i]).select_cell(
                        pos, dds) == 1:
                    return 1
        return 0

    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef int select_point(self, np.float64_t pos[3]) nogil:
        with gil:
            for i in range(self.n_obj):
                if (<SelectorObject>self.selectors[i]).select_point(pos) == 1:
                    return 1
        return 0

    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef int select_sphere(self, np.float64_t pos[3], np.float64_t radius) nogil:
        with gil:
            for i in range(self.n_obj):
                if (<SelectorObject>self.selectors[i]).select_sphere(
                        pos, radius) == 1:
                    return 1
        return 0

    def _hash_vals(self):
        v = ("chained_or",)
        for s in self.selectors:
            v += s._hash_vals()
        return v

union_selector = ChainedBooleanORSelector
