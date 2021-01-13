cdef class PointSelector(SelectorObject):
    cdef public np.float64_t p[3]

    def __init__(self, dobj):
        cdef np.float64_t[:] DLE = _ensure_code(dobj.ds.domain_left_edge)
        cdef np.float64_t[:] DRE = _ensure_code(dobj.ds.domain_right_edge)
        for i in range(3):
            self.p[i] = _ensure_code(dobj.p[i])

            # ensure the point lies in the domain
            if self.periodicity[i]:
                self.p[i] = np.fmod(self.p[i], self.domain_width[i])
                if self.p[i] < DLE[i]:
                    self.p[i] += self.domain_width[i]
                elif self.p[i] >= DRE[i]:
                    self.p[i] -= self.domain_width[i]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        if (pos[0] - 0.5*dds[0] <= self.p[0] < pos[0]+0.5*dds[0] and
            pos[1] - 0.5*dds[1] <= self.p[1] < pos[1]+0.5*dds[1] and
            pos[2] - 0.5*dds[2] <= self.p[2] < pos[2]+0.5*dds[2]):
            return 1
        else:
            return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_sphere(self, np.float64_t pos[3], np.float64_t radius) nogil:
        cdef int i
        cdef np.float64_t dist, dist2 = 0
        for i in range(3):
            dist = self.periodic_difference(pos[i], self.p[i], i)
            dist2 += dist*dist
        if dist2 <= radius*radius: return 1
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_bbox(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        # point definitely can only be in one cell
        if (left_edge[0] <= self.p[0] < right_edge[0] and
            left_edge[1] <= self.p[1] < right_edge[1] and
            left_edge[2] <= self.p[2] < right_edge[2]):
            return 1
        else:
            return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_bbox_edge(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        # point definitely can only be in one cell
        # Return 2 in all cases to indicate that the point only overlaps
        # portion of box
        if (left_edge[0] <= self.p[0] <= right_edge[0] and
            left_edge[1] <= self.p[1] <= right_edge[1] and
            left_edge[2] <= self.p[2] <= right_edge[2]):
            return 2
        else:
            return 0

    def _hash_vals(self):
        return (("p[0]", self.p[0]),
                ("p[1]", self.p[1]),
                ("p[2]", self.p[2]))

    def _get_state_attnames(self):
        return ('p', )

point_selector = PointSelector
