from libc.math cimport sin, cos, atan2, sqrt, acos

cdef class CuttingPlaneSelector(SelectorObject):
    cdef public np.float64_t norm_vec[3]
    cdef public np.float64_t d

    def __init__(self, dobj):
        cdef int i
        for i in range(3):
            self.norm_vec[i] = dobj._norm_vec[i]
        self.d = _ensure_code(dobj._d)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) noexcept nogil:
        cdef np.float64_t left_edge[3]
        cdef np.float64_t right_edge[3]
        cdef int i
        for i in range(3):
            left_edge[i] = pos[i] - 0.5*dds[i]
            right_edge[i] = pos[i] + 0.5*dds[i]
        return self.select_bbox(left_edge, right_edge)

    cdef int select_point(self, np.float64_t pos[3]) noexcept nogil:
        # two 0-volume constructs don't intersect
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_sphere(self, np.float64_t pos[3], np.float64_t radius) noexcept nogil:
        cdef int i
        cdef np.float64_t height = self.d
        for i in range(3) :
            height += pos[i] * self.norm_vec[i]
        if height*height <= radius*radius : return 1
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_bbox(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) noexcept nogil:
        cdef int i, j, k, n
        cdef np.float64_t *arr[2]
        cdef np.float64_t pos[3]
        cdef np.float64_t gd
        arr[0] = left_edge
        arr[1] = right_edge
        all_under = 1
        all_over = 1
        # Check each corner
        for i in range(2):
            pos[0] = arr[i][0]
            for j in range(2):
                pos[1] = arr[j][1]
                for k in range(2):
                    pos[2] = arr[k][2]
                    gd = self.d
                    for n in range(3):
                        gd += pos[n] * self.norm_vec[n]
                    # this allows corners and faces on the low-end to
                    # collide, while not selecting cells on the high-side
                    if i == 0 and j == 0 and k == 0 :
                        if gd <= 0: all_over = 0
                        if gd >= 0: all_under = 0
                    else :
                        if gd < 0: all_over = 0
                        if gd > 0: all_under = 0
        if all_over == 1 or all_under == 1:
            return 0
        return 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_bbox_edge(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) noexcept nogil:
        cdef int i, j, k, n
        cdef np.float64_t *arr[2]
        cdef np.float64_t pos[3]
        cdef np.float64_t gd
        arr[0] = left_edge
        arr[1] = right_edge
        all_under = 1
        all_over = 1
        # Check each corner
        for i in range(2):
            pos[0] = arr[i][0]
            for j in range(2):
                pos[1] = arr[j][1]
                for k in range(2):
                    pos[2] = arr[k][2]
                    gd = self.d
                    for n in range(3):
                        gd += pos[n] * self.norm_vec[n]
                    # this allows corners and faces on the low-end to
                    # collide, while not selecting cells on the high-side
                    if i == 0 and j == 0 and k == 0 :
                        if gd <= 0: all_over = 0
                        if gd >= 0: all_under = 0
                    else :
                        if gd < 0: all_over = 0
                        if gd > 0: all_under = 0
        if all_over == 1 or all_under == 1:
            return 0
        return 2 # a box of non-zeros volume can't be inside a plane

    def _hash_vals(self):
        return (("norm_vec[0]", self.norm_vec[0]),
                ("norm_vec[1]", self.norm_vec[1]),
                ("norm_vec[2]", self.norm_vec[2]),
                ("d", self.d))

    def _get_state_attnames(self):
        return ("d", "norm_vec")

cdef class SphericalCuttingPlaneSelector(CuttingPlaneSelector):

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void transform_rtp_to_xyz(self, np.float64_t in_pos[3], np.float64_t out_pos[3]) noexcept nogil:
        # See `spherical_coordinates.py` for more details
        # r => in_pos[0] theta => in_pos[1] phi => in_pos[2]
        out_pos[0] = in_pos[0] * cos(in_pos[2]) * sin(in_pos[1])
        out_pos[1] = in_pos[0] * sin(in_pos[2]) * sin(in_pos[1])
        out_pos[2] = in_pos[0] * cos(in_pos[1])

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void transform_xyz_to_rtp(self, np.float64_t in_pos[3], np.float64_t out_pos[3]) noexcept nogil:
        # we have to have 012 be xyz and 012 be rtp
        out_pos[0] = sqrt(in_pos[0]*in_pos[0] + in_pos[1]*in_pos[1] + in_pos[2]*in_pos[2])
        out_pos[1] = acos(in_pos[2] / out_pos[0])
        out_pos[2] = atan2(in_pos[1], in_pos[0])

cutting_selector = CuttingPlaneSelector

spherical_cutting_selector = SphericalCuttingPlaneSelector
