from yt.utilities.lib.coordinate_utilities cimport spherical_to_cartesian, cartesian_to_spherical
from numpy.math cimport PI as NPY_PI

cdef class CuttingPlaneSelector(SelectorObject):
    cdef public np.float64_t norm_vec[3]  # the unit-normal for the plane
    cdef public np.float64_t d  # the shortest distance from plane to origin

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
        # the bbox selection here works by calculating the signed-distance from
        # the plane to each vertex of the bounding box. If there is no
        # intersection, the signed-distance for every vertex will have the same
        # sign whereas if the sign flips then the plane must intersect the
        # bounding box.
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


cdef class CuttingPlaneTransformed(CuttingPlaneSelector):
    # a base class for cartesian cutting planes through data that is not
    # in cartesian coordinates.

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void transform_vertex_pos(self, np.float64_t pos_in[3], np.float64_t pos_out[3]) noexcept nogil:
        # child class must implement: must transform from dataset native
        # coordinates to cartesian coordinates (with (x,y,z) ordering)
        pass

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_bbox(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) noexcept nogil:
        # child classes may over-ride if needed
        cdef np.int64_t n_points[3]
        n_points[0] = 2
        n_points[1] = 2
        n_points[2] = 2
        return self._select_bbox(left_edge, right_edge, n_points)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int _select_bbox(self,
                          np.float64_t left_edge[3],
                          np.float64_t right_edge[3],
                          np.int64_t n_points[3],
                          ) noexcept nogil:

        # the bbox selection here works by calculating the signed-distance from
        # the plane to each vertex of the bounding box. If there is no
        # intersection, the signed-distance for every vertex will have the same
        # sign whereas if the sign flips then the plane must intersect the
        # bounding box.
        #
        # left_edge, right_edge
        #   the left/right bounds of the bounding box
        # n_points :
        #   the number of points to check in each dimension.
        #   (2,2,2) is equivalent to checking the corners of the
        #   bounding box.

        cdef int i, j, k, n
        cdef np.float64_t pos[3]
        cdef np.float64_t dpos[3]
        cdef np.float64_t pos_cart[3]
        cdef np.float64_t gd, n_pts

        for i in range(3):
            dpos[i] = (right_edge[i] - left_edge[i]) / (<float> n_points[i] - 1.0)

        all_under = 1
        all_over = 1

        for i in range(n_points[0]):
            pos[0] = left_edge[0] + <float> i * dpos[0]
            for j in range(n_points[1]):
                pos[1] = left_edge[1] + <float> j * dpos[1]
                for k in range(n_points[2]):
                    pos[2] = left_edge[2] + <float> k * dpos[2]
                    self.transform_vertex_pos(pos, pos_cart)
                    gd = self.d
                    for n in range(3):
                        gd += pos_cart[n] * self.norm_vec[n]
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


cdef class SphericalCuttingPlaneSelector(CuttingPlaneTransformed):

    # intersection of a cartesian plane with data in spherical coordinates.
    # expected ordering is (r, theta, phi), where theta is the colatitude
    # angle (bounds of 0 to pi) and phi is the azimuthal/longitudinal
    # angle (bounds 0 to 2pi).

    cdef public np.float64_t r_min # the minimum radius for possible intersection
    cdef public np.float64_t c_rtp[3]

    def __init__(self, dobj):

        cdef np.float64_t xyz[3]
        cdef int i
        super().__init__(dobj)

        # any points at r < |d|, where d is the minimum distance-vector to the
        # plane, cannot intersect the plane. Record r_min here for convenience:
        self.r_min = fabs(self.d)

        # also record the spherical coordinates of the point on the plane
        # closest to the origin
        for i in range(3):
            xyz[i] = - self.norm_vec[i] * self.d  # cartesian position
        self.c_rtp[0],  self.c_rtp[1],  self.c_rtp[2] = cartesian_to_spherical(xyz[0], xyz[1], xyz[2])

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void transform_vertex_pos(self, np.float64_t pos_in[3], np.float64_t pos_out[3]) noexcept nogil:
        # r => in_pos[0] theta => in_pos[1] phi => in_pos[2]
        pos_out[0], pos_out[1], pos_out[2] = spherical_to_cartesian(pos_in[0], pos_in[1], pos_in[2])


    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_bbox(self,
                         np.float64_t left_edge[3],
                         np.float64_t right_edge[3]) noexcept nogil:
         # left/right edge here are in spherical coordinates in (r, theta, phi)

         cdef int selected, idim
         cdef np.float64_t left_edge_c[3], right_edge_c[3],
         cdef np.int64_t n_points[3]
         cdef np.float64_t NPY_PI_4 = NPY_PI / 4.0
         cdef np.float64_t dangle, PI2


         # set the number of points to check depending on angular range
         # of element: small elements will only check the bounding box
         # verts as normal. elements with large angular range, > pi/4, will
         # add additional vertices.
         for idim in range(3):
            n_points[idim] = 2
         # check the angular ranges and add extra points if too large...
         for idim in range(1,3):
            dangle = right_edge[idim] - left_edge[idim]
            if dangle >= NPY_PI_4:
                n_points[idim] = 2 + <int> (dangle / NPY_PI_4)

         # run the plane-vertex distance check (vertex positions are converted to
         # cartesian within _select_bbox).
         selected = self._select_bbox(left_edge, right_edge, n_points)

         if selected == 0:
            # there is one more special case to consider!
            # When all vertices lie on one side of the plane, intersection
            # is still possible if the plane intersects the outer cusp of the
            # spherical volume element. **BUT** if we've reached this far,
            # the only way for this to happen is if the position of the point
            # on the plane that is closest to the origin lies within the
            # element itself.
            if self.c_rtp[0] > left_edge[0]:
                for idim in range(1,3):
                    if self.c_rtp[idim] <= left_edge[idim]:
                        return 0
                    if self.c_rtp[idim] >= right_edge[idim]:
                        return 0
                return 1

         return selected

    def _select_single_bbox(self,
                  left_edge_in,
                  right_edge_in):

         # useful for direct testing without having to initialize
         # full yt data objects

         cdef np.float64_t left_edge[3]
         cdef np.float64_t right_edge[3]

         for i in range(3):
             left_edge[i] = left_edge_in[i]
             right_edge[i] = right_edge_in[i]

         return self.select_bbox(left_edge, right_edge)


cutting_selector = CuttingPlaneSelector

cutting_mixed_spherical_selector = SphericalCuttingPlaneSelector
