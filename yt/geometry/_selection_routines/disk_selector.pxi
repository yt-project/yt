cdef class DiskSelector(SelectorObject):
    cdef public np.float64_t norm_vec[3]
    cdef public np.float64_t center[3]
    cdef public np.float64_t radius, radius2
    cdef public np.float64_t height

    def __init__(self, dobj):
        cdef int i
        for i in range(3):
            self.norm_vec[i] = dobj._norm_vec[i]
            self.center[i] = _ensure_code(dobj.center[i])
        self.radius = _ensure_code(dobj.radius)
        self.radius2 = self.radius * self.radius
        self.height = _ensure_code(dobj.height)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        return self.select_point(pos)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_point(self, np.float64_t pos[3]) nogil:
        cdef np.float64_t h, d, r2, temp
        cdef int i
        h = d = 0
        for i in range(3):
            temp = self.periodic_difference(pos[i], self.center[i], i)
            h += temp * self.norm_vec[i]
            d += temp*temp
        r2 = (d - h*h)
        if fabs(h) <= self.height and r2 <= self.radius2: return 1
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_sphere(self, np.float64_t pos[3], np.float64_t radius) nogil:
        cdef np.float64_t h, d, r2, temp
        cdef int i
        h = d = 0
        for i in range(3):
            temp = self.periodic_difference(pos[i], self.center[i], i)
            h += temp * self.norm_vec[i]
            d += temp*temp
        r2 = (d - h*h)
        d = self.radius+radius
        if fabs(h) <= self.height+radius and r2 <= d*d: return 1
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_bbox(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        # Until we can get our OBB/OBB intersection correct, disable this.
        return 1
        # cdef np.float64_t *arr[2]
        # cdef np.float64_t pos[3]
        # cdef np.float64_t H, D, R2, temp
        # cdef int i, j, k, n
        # cdef int all_under = 1
        # cdef int all_over = 1
        # cdef int any_radius = 0
        # # A moment of explanation (revised):
        # #    The disk and bounding box collide if any of the following are true:
        # #    1) the center of the disk is inside the bounding box
        # #    2) any corner of the box lies inside the disk
        # #    3) the box spans the plane (!all_under and !all_over) and at least
        # #       one corner is within the cylindrical radius

        # # check if disk center lies inside bbox
        # if left_edge[0] <= self.center[0] <= right_edge[0] and \
        #    left_edge[1] <= self.center[1] <= right_edge[1] and \
        #    left_edge[2] <= self.center[2] <= right_edge[2] :
        #     return 1

        # # check all corners
        # arr[0] = left_edge
        # arr[1] = right_edge
        # for i in range(2):
        #     pos[0] = arr[i][0]
        #     for j in range(2):
        #         pos[1] = arr[j][1]
        #         for k in range(2):
        #             pos[2] = arr[k][2]
        #             H = D = 0
        #             for n in range(3):
        #                 temp = self.difference(pos[n], self.center[n], n)
        #                 H += (temp * self.norm_vec[n])
        #                 D += temp*temp
        #             R2 = (D - H*H)
        #             if R2 < self.radius2 :
        #                 any_radius = 1
        #                 if fabs(H) < self.height: return 1
        #             if H < 0: all_over = 0
        #             if H > 0: all_under = 0
        # if all_over == 0 and all_under == 0 and any_radius == 1: return 1
        # return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_bbox_edge(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        # Until we can get our OBB/OBB intersection correct, disable this.
        return 2
        # cdef np.float64_t *arr[2]
        # cdef np.float64_t pos[3], H, D, R2, temp
        # cdef int i, j, k, n
        # cdef int all_under = 1
        # cdef int all_over = 1
        # cdef int any_radius = 0
        # # A moment of explanation (revised):
        # #    The disk and bounding box collide if any of the following are true:
        # #    1) the center of the disk is inside the bounding box
        # #    2) any corner of the box lies inside the disk
        # #    3) the box spans the plane (!all_under and !all_over) and at least
        # #       one corner is within the cylindrical radius

        # # check if disk center lies inside bbox
        # if left_edge[0] <= self.center[0] <= right_edge[0] and \
        #    left_edge[1] <= self.center[1] <= right_edge[1] and \
        #    left_edge[2] <= self.center[2] <= right_edge[2] :
        #     return 1

        # # check all corners
        # arr[0] = left_edge
        # arr[1] = right_edge
        # for i in range(2):
        #     pos[0] = arr[i][0]
        #     for j in range(2):
        #         pos[1] = arr[j][1]
        #         for k in range(2):
        #             pos[2] = arr[k][2]
        #             H = D = 0
        #             for n in range(3):
        #                 temp = self.periodic_difference(
        #                     pos[n], self.center[n], n)
        #                 H += (temp * self.norm_vec[n])
        #                 D += temp*temp
        #             R2 = (D - H*H)
        #             if R2 < self.radius2 :
        #                 any_radius = 1
        #                 if fabs(H) < self.height: return 1
        #             if H < 0: all_over = 0
        #             if H > 0: all_under = 0
        # if all_over == 0 and all_under == 0 and any_radius == 1: return 1
        # return 0

    def _hash_vals(self):
        return (("norm_vec[0]", self.norm_vec[0]),
                ("norm_vec[1]", self.norm_vec[1]),
                ("norm_vec[2]", self.norm_vec[2]),
                ("center[0]", self.center[0]),
                ("center[1]", self.center[1]),
                ("center[2]", self.center[2]),
                ("radius", self.radius),
                ("radius2", self.radius2),
                ("height", self.height))

    def _get_state_attnames(self):
        return ("radius", "radius2", "height", "norm_vec", "center")


disk_selector = DiskSelector
