cdef class EllipsoidSelector(SelectorObject):
    cdef public np.float64_t vec[3][3]
    cdef public np.float64_t mag[3]
    cdef public np.float64_t center[3]

    def __init__(self, dobj):
        cdef int i
        _ensure_code(dobj.center)
        _ensure_code(dobj._e0)
        _ensure_code(dobj._e1)
        _ensure_code(dobj._e2)
        _ensure_code(dobj._A)
        _ensure_code(dobj._B)
        _ensure_code(dobj._C)
        for i in range(3):
            self.center[i] = dobj.center[i]
            self.vec[0][i] = dobj._e0[i]
            self.vec[1][i] = dobj._e1[i]
            self.vec[2][i] = dobj._e2[i]
        self.mag[0] = dobj._A
        self.mag[1] = dobj._B
        self.mag[2] = dobj._C

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        return self.select_point(pos)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_point(self, np.float64_t pos[3]) nogil:
        cdef np.float64_t dot_evec[3]
        cdef np.float64_t dist
        cdef int i, j
        dot_evec[0] = dot_evec[1] = dot_evec[2] = 0
        # Calculate the rotated dot product
        for i in range(3): # axis
            dist = self.periodic_difference(pos[i], self.center[i], i)
            for j in range(3):
                dot_evec[j] += dist * self.vec[j][i]
        dist = 0.0
        for i in range(3):
            dist += (dot_evec[i] * dot_evec[i])/(self.mag[i] * self.mag[i])
        if dist <= 1.0: return 1
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_sphere(self, np.float64_t pos[3], np.float64_t radius) nogil:
        # this is the sphere selection
        cdef int i
        cdef np.float64_t dist, dist2_max, dist2 = 0
        for i in range(3):
            dist = self.periodic_difference(pos[i], self.center[i], i)
            dist2 += dist * dist
        dist2_max = (self.mag[0] + radius) * (self.mag[0] + radius)
        if dist2 <= dist2_max:
            return 1
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_bbox(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        # This is the sphere selection
        cdef int i
        cdef np.float64_t box_center, relcenter, closest, dist, edge, dist_max
        if left_edge[0] <= self.center[0] <= right_edge[0] and \
           left_edge[1] <= self.center[1] <= right_edge[1] and \
           left_edge[2] <= self.center[2] <= right_edge[2]:
            return 1
        # http://www.gamedev.net/topic/335465-is-this-the-simplest-sphere-aabb-collision-test/
        dist = 0
        for i in range(3):
            box_center = (right_edge[i] + left_edge[i])/2.0
            relcenter = self.periodic_difference(box_center, self.center[i], i)
            edge = right_edge[i] - left_edge[i]
            closest = relcenter - fclip(relcenter, -edge/2.0, edge/2.0)
            dist += closest * closest
        dist_max = self.mag[0] * self.mag[0]
        if dist <= dist_max:
            return 1
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_bbox_edge(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        # This is the sphere selection
        cdef int i
        cdef np.float64_t box_center, relcenter, closest, farthest, cdist, fdist, edge
        if left_edge[0] <= self.center[0] <= right_edge[0] and \
           left_edge[1] <= self.center[1] <= right_edge[1] and \
           left_edge[2] <= self.center[2] <= right_edge[2]:
            fdist = 0
            for i in range(3):
                edge = right_edge[i] - left_edge[i]
                box_center = (right_edge[i] + left_edge[i])/2.0
                relcenter = self.periodic_difference(
                    box_center, self.center[i], i)
                farthest = relcenter + fclip(relcenter, -edge/2.0, edge/2.0)
                fdist += farthest*farthest
                if fdist >= self.mag[0]**2: return 2
            return 1
        # http://www.gamedev.net/topic/335465-is-this-the-simplest-sphere-aabb-collision-test/
        cdist = 0
        fdist = 0
        for i in range(3):
            box_center = (right_edge[i] + left_edge[i])/2.0
            relcenter = self.periodic_difference(box_center, self.center[i], i)
            edge = right_edge[i] - left_edge[i]
            closest = relcenter - fclip(relcenter, -edge/2.0, edge/2.0)
            farthest = relcenter + fclip(relcenter, -edge/2.0, edge/2.0)
            cdist += closest * closest
            fdist += farthest * farthest
            if cdist > self.mag[0]**2: return 0
        if fdist < self.mag[0]**2:
            return 1
        else:
            return 2

    def _hash_vals(self):
        return (("vec[0][0]", self.vec[0][0]),
                ("vec[0][1]", self.vec[0][1]),
                ("vec[0][2]", self.vec[0][2]),
                ("vec[1][0]", self.vec[1][0]),
                ("vec[1][1]", self.vec[1][1]),
                ("vec[1][2]", self.vec[1][2]),
                ("vec[2][0]", self.vec[2][0]),
                ("vec[2][1]", self.vec[2][1]),
                ("vec[2][2]", self.vec[2][2]),
                ("mag[0]", self.mag[0]),
                ("mag[1]", self.mag[1]),
                ("mag[2]", self.mag[2]),
                ("center[0]", self.center[0]),
                ("center[1]", self.center[1]),
                ("center[2]", self.center[2]))

    def _get_state_attnames(self):
        return ("mag", "center", "vec")


ellipsoid_selector = EllipsoidSelector
