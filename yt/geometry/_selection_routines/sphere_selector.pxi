cdef class SphereSelector(SelectorObject):
    cdef public np.float64_t radius
    cdef public np.float64_t radius2
    cdef public np.float64_t center[3]
    cdef np.float64_t bbox[3][2]
    cdef public bint check_box[3]

    def __init__(self, dobj):
        for i in range(3):
            self.center[i] = _ensure_code(dobj.center[i])
        self.radius = _ensure_code(dobj.radius)
        self.radius2 = self.radius * self.radius
        self.set_bbox(_ensure_code(dobj.center))
        for i in range(3):
            if self.bbox[i][0] < dobj.ds.domain_left_edge[i]:
                self.check_box[i] = False
            elif self.bbox[i][1] > dobj.ds.domain_right_edge[i]:
                self.check_box[i] = False
            else:
                self.check_box[i] = True

    def set_bbox(self, center):
        for i in range(3):
            self.center[i] = center[i]
            self.bbox[i][0] = self.center[i] - self.radius
            self.bbox[i][1] = self.center[i] + self.radius

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        # sphere center either inside cell or center of cell lies inside sphere
        if (pos[0] - 0.5*dds[0] <= self.center[0] <= pos[0]+0.5*dds[0] and
            pos[1] - 0.5*dds[1] <= self.center[1] <= pos[1]+0.5*dds[1] and
            pos[2] - 0.5*dds[2] <= self.center[2] <= pos[2]+0.5*dds[2]):
            return 1
        return self.select_point(pos)
        # # langmm: added to allow sphere to intersect edge/corner of cell
        # cdef np.float64_t LE[3]
        # cdef np.float64_t RE[3]
        # cdef int i
        # for i in range(3):
        #     LE[i] = pos[i] - 0.5*dds[i]
        #     RE[i] = pos[i] + 0.5*dds[i]
        # return self.select_bbox(LE, RE)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_point(self, np.float64_t pos[3]) nogil:
        cdef int i
        cdef np.float64_t dist, dist2 = 0
        for i in range(3):
            if self.check_box[i] and \
              (pos[i] < self.bbox[i][0] or
               pos[i] > self.bbox[i][1]):
                return 0
            dist = _periodic_dist(pos[i], self.center[i], self.domain_width[i],
                                  self.periodicity[i])
            dist2 += dist*dist
            if dist2 > self.radius2: return 0
        return 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_sphere(self, np.float64_t pos[3], np.float64_t radius) nogil:
        cdef int i
        cdef np.float64_t dist, dist2 = 0
        for i in range(3):
            dist = self.periodic_difference(pos[i], self.center[i], i)
            dist2 += dist*dist
        dist = self.radius+radius
        if dist2 <= dist*dist: return 1
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_bbox(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        cdef np.float64_t box_center, relcenter, closest, dist, edge
        cdef int i
        if (left_edge[0] <= self.center[0] < right_edge[0] and
            left_edge[1] <= self.center[1] < right_edge[1] and
            left_edge[2] <= self.center[2] < right_edge[2]):
            return 1
        for i in range(3):
            if not self.check_box[i]: continue
            if right_edge[i] < self.bbox[i][0] or \
               left_edge[i] > self.bbox[i][1]:
                return 0
        # http://www.gamedev.net/topic/335465-is-this-the-simplest-sphere-aabb-collision-test/
        dist = 0
        for i in range(3):
            # Early terminate
            box_center = (right_edge[i] + left_edge[i])/2.0
            relcenter = self.periodic_difference(box_center, self.center[i], i)
            edge = right_edge[i] - left_edge[i]
            closest = relcenter - fclip(relcenter, -edge/2.0, edge/2.0)
            dist += closest*closest
            if dist > self.radius2: return 0
        return 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_bbox_edge(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        cdef np.float64_t box_center, relcenter, closest, farthest, cdist, fdist, edge
        cdef int i
        if (left_edge[0] <= self.center[0] <= right_edge[0] and
            left_edge[1] <= self.center[1] <= right_edge[1] and
            left_edge[2] <= self.center[2] <= right_edge[2]):
            fdist = 0
            for i in range(3):
                edge = right_edge[i] - left_edge[i]
                box_center = (right_edge[i] + left_edge[i])/2.0
                relcenter = self.periodic_difference(
                    box_center, self.center[i], i)
                if relcenter >= 0:
                    farthest = relcenter + edge/2.0
                else:
                    farthest = relcenter - edge/2.0
                # farthest = relcenter + fclip(relcenter, -edge/2.0, edge/2.0)
                fdist += farthest*farthest
                if fdist >= self.radius2:
                    return 2  # Box extends outside sphere
            return 1  # Box entirely inside sphere
        for i in range(3):
            if not self.check_box[i]: continue
            if right_edge[i] < self.bbox[i][0] or \
               left_edge[i] > self.bbox[i][1]:
                return 0  # Box outside sphere bounding box
        # http://www.gamedev.net/topic/335465-is-this-the-simplest-sphere-aabb-collision-test/
        cdist = 0
        fdist = 0
        for i in range(3):
            # Early terminate
            box_center = (right_edge[i] + left_edge[i])/2.0
            relcenter = self.periodic_difference(box_center, self.center[i], i)
            edge = right_edge[i] - left_edge[i]
            closest = relcenter - fclip(relcenter, -edge/2.0, edge/2.0)
            if relcenter >= 0:
                farthest = relcenter + edge/2.0
            else:
                farthest = relcenter - edge/2.0
            #farthest = relcenter + fclip(relcenter, -edge/2.0, edge/2.0)
            cdist += closest*closest
            fdist += farthest*farthest
            if cdist > self.radius2:
                return 0  # Box does not overlap sphere
        if fdist < self.radius2:
            return 1  # Sphere extends to entirely contain box
        else:
            return 2  # Sphere only partially overlaps box

    def _hash_vals(self):
        return (("radius", self.radius),
                ("radius2", self.radius2),
                ("center[0]", self.center[0]),
                ("center[1]", self.center[1]),
                ("center[2]", self.center[2]))

    def _get_state_attnames(self):
        return ("radius", "radius2", "center", "check_box")

    def __setstate__(self, hashes):
        super(SphereSelector, self).__setstate__(hashes)
        self.set_bbox(self.center)


sphere_selector = SphereSelector
