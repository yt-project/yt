cdef struct IntegrationAccumulator:
    np.float64_t *t
    np.float64_t *dt
    np.uint8_t *child_mask
    int hits

cdef void dt_sampler(
             VolumeContainer *vc,
             np.float64_t v_pos[3],
             np.float64_t v_dir[3],
             np.float64_t enter_t,
             np.float64_t exit_t,
             int index[3],
             void *data) nogil:
    cdef IntegrationAccumulator *am = <IntegrationAccumulator *> data
    cdef int di = (index[0]*vc.dims[1]+index[1])*vc.dims[2]+index[2]
    if am.child_mask[di] == 0 or enter_t == exit_t:
        return
    am.hits += 1
    am.t[di] = enter_t
    am.dt[di] = (exit_t - enter_t)

cdef class RaySelector(SelectorObject):

    cdef public np.float64_t p1[3]
    cdef public np.float64_t p2[3]
    cdef public np.float64_t vec[3]

    def __init__(self, dobj):
        cdef int i
        _ensure_code(dobj.start_point)
        _ensure_code(dobj.end_point)
        for i in range(3):
            self.vec[i] = dobj.vec[i]
            self.p1[i] = dobj.start_point[i]
            self.p2[i] = dobj.end_point[i]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def fill_mask(self, gobj):
        cdef np.ndarray[np.float64_t, ndim=3] t, dt
        cdef np.ndarray[np.uint8_t, ndim=3, cast=True] child_mask
        cdef int i
        cdef int total = 0
        cdef IntegrationAccumulator *ia
        ia = <IntegrationAccumulator *> malloc(sizeof(IntegrationAccumulator))
        cdef VolumeContainer vc
        mask = np.zeros(gobj.ActiveDimensions, dtype='uint8')
        t = np.zeros(gobj.ActiveDimensions, dtype="float64")
        dt = np.zeros(gobj.ActiveDimensions, dtype="float64") - 1
        child_mask = gobj.child_mask
        ia.t = <np.float64_t *> t.data
        ia.dt = <np.float64_t *> dt.data
        ia.child_mask = <np.uint8_t *> child_mask.data
        ia.hits = 0
        _ensure_code(gobj.LeftEdge)
        _ensure_code(gobj.RightEdge)
        _ensure_code(gobj.dds)
        for i in range(3):
            vc.left_edge[i] = gobj.LeftEdge[i]
            vc.right_edge[i] = gobj.RightEdge[i]
            vc.dds[i] = gobj.dds[i]
            vc.idds[i] = 1.0/gobj.dds[i]
            vc.dims[i] = dt.shape[i]
        walk_volume(&vc, self.p1, self.vec, dt_sampler, <void*> ia)
        for i in range(dt.shape[0]):
            for j in range(dt.shape[1]):
                for k in range(dt.shape[2]):
                    if dt[i, j, k] >= 0:
                        mask[i, j, k] = 1
                        total += 1
        free(ia)
        if total == 0: return None
        return mask.astype("bool")

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def get_dt(self, gobj):
        cdef np.ndarray[np.float64_t, ndim=3] t, dt
        cdef np.ndarray[np.float64_t, ndim=1] tr, dtr
        cdef np.ndarray[np.uint8_t, ndim=3, cast=True] child_mask
        cdef int i, j, k, ni
        cdef IntegrationAccumulator *ia
        ia = <IntegrationAccumulator *> malloc(sizeof(IntegrationAccumulator))
        cdef VolumeContainer vc
        t = np.zeros(gobj.ActiveDimensions, dtype="float64")
        dt = np.zeros(gobj.ActiveDimensions, dtype="float64") - 1
        child_mask = gobj.child_mask
        ia.t = <np.float64_t *> t.data
        ia.dt = <np.float64_t *> dt.data
        ia.child_mask = <np.uint8_t *> child_mask.data
        ia.hits = 0
        _ensure_code(gobj.LeftEdge)
        _ensure_code(gobj.RightEdge)
        _ensure_code(gobj.dds)
        for i in range(3):
            vc.left_edge[i] = gobj.LeftEdge[i]
            vc.right_edge[i] = gobj.RightEdge[i]
            vc.dds[i] = gobj.dds[i]
            vc.idds[i] = 1.0/gobj.dds[i]
            vc.dims[i] = dt.shape[i]
        walk_volume(&vc, self.p1, self.vec, dt_sampler, <void*> ia)
        tr = np.zeros(ia.hits, dtype="float64")
        dtr = np.zeros(ia.hits, dtype="float64")
        ni = 0
        for i in range(dt.shape[0]):
            for j in range(dt.shape[1]):
                for k in range(dt.shape[2]):
                    if dt[i, j, k] >= 0:
                        tr[ni] = t[i, j, k]
                        dtr[ni] = dt[i, j, k]
                        ni += 1
        if not (ni == ia.hits):
            print(ni, ia.hits)
        free(ia)
        return dtr, tr

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def get_dt_mesh(self, mesh, nz, int offset):
        cdef np.ndarray[np.float64_t, ndim=3] t, dt
        cdef np.ndarray[np.float64_t, ndim=1] tr, dtr
        cdef np.ndarray[np.uint8_t, ndim=3, cast=True] child_mask
        cdef int i, j, k, ni
        cdef np.float64_t LE[3]
        cdef np.float64_t RE[3]
        cdef np.float64_t pos
        cdef IntegrationAccumulator *ia
        ia = <IntegrationAccumulator *> malloc(sizeof(IntegrationAccumulator))
        cdef np.ndarray[np.float64_t, ndim=2] coords
        cdef np.ndarray[np.int64_t, ndim=2] indices
        indices = mesh.connectivity_indices
        coords = _ensure_code(mesh.connectivity_coords)
        cdef int nc = indices.shape[0]
        cdef int nv = indices.shape[1]
        if nv != 8:
            raise NotImplementedError
        cdef VolumeContainer vc
        child_mask = np.ones((1,1,1), dtype="uint8")
        t = np.zeros((1,1,1), dtype="float64")
        dt = np.zeros((1,1,1), dtype="float64") - 1
        tr = np.zeros(nz, dtype="float64")
        dtr = np.zeros(nz, dtype="float64")
        ia.t = <np.float64_t *> t.data
        ia.dt = <np.float64_t *> dt.data
        ia.child_mask = <np.uint8_t *> child_mask.data
        ia.hits = 0
        ni = 0
        for i in range(nc):
            for j in range(3):
                LE[j] = 1e60
                RE[j] = -1e60
            for j in range(nv):
                for k in range(3):
                    pos = coords[indices[i, j] - offset, k]
                    LE[k] = fmin(pos, LE[k])
                    RE[k] = fmax(pos, RE[k])
            for j in range(3):
                vc.left_edge[j] = LE[j]
                vc.right_edge[j] = RE[j]
                vc.dds[j] = RE[j] - LE[j]
                vc.idds[j] = 1.0/vc.dds[j]
                vc.dims[j] = 1
            t[0,0,0] = dt[0,0,0] = -1
            walk_volume(&vc, self.p1, self.vec, dt_sampler, <void*> ia)
            if dt[0,0,0] >= 0:
                tr[ni] = t[0,0,0]
                dtr[ni] = dt[0,0,0]
                ni += 1
        free(ia)
        return dtr, tr

    cdef int select_point(self, np.float64_t pos[3]) nogil:
        # two 0-volume constructs don't intersect
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_sphere(self, np.float64_t pos[3], np.float64_t radius) nogil:

        cdef int i
        cdef np.float64_t length = norm(self.vec)
        cdef np.float64_t r[3]
        for i in range(3):
            r[i] = pos[i] - self.p1[i]
        # the projected position of the sphere along the ray
        cdef np.float64_t l = dot(r, self.vec) / length
        # the square of the impact parameter
        cdef np.float64_t b_sqr = dot(r, r) - l*l

        # only accept spheres with radii larger than the impact parameter and
        # with a projected position along the ray no more than a radius away
        # from the ray
        if -radius < l and l < (length+radius) and b_sqr < radius*radius:
            return 1

        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_bbox(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        cdef int i, rv
        cdef VolumeContainer vc
        cdef IntegrationAccumulator *ia
        ia = <IntegrationAccumulator *> malloc(sizeof(IntegrationAccumulator))
        cdef np.float64_t dt[1]
        cdef np.float64_t t[1]
        cdef np.uint8_t cm[1]
        for i in range(3):
            vc.left_edge[i] = left_edge[i]
            vc.right_edge[i] = right_edge[i]
            vc.dds[i] = right_edge[i] - left_edge[i]
            vc.idds[i] = 1.0/vc.dds[i]
            vc.dims[i] = 1
        t[0] = dt[0] = 0.0
        cm[0] = 1
        ia.t = t
        ia.dt = dt
        ia.child_mask = cm
        ia.hits = 0
        walk_volume(&vc, self.p1, self.vec, dt_sampler, <void*> ia)
        rv = 0
        if ia.hits > 0:
            rv = 1
        free(ia)
        return rv

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_bbox_edge(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        cdef int i
        cdef np.uint8_t cm = 1
        cdef VolumeContainer vc
        cdef IntegrationAccumulator ia
        cdef np.float64_t dt, t
        for i in range(3):
            vc.left_edge[i] = left_edge[i]
            vc.right_edge[i] = right_edge[i]
            vc.dds[i] = right_edge[i] - left_edge[i]
            vc.idds[i] = 1.0/vc.dds[i]
            vc.dims[i] = 1
        t = dt = 0.0
        ia.t = &t
        ia.dt = &dt
        ia.child_mask = &cm
        ia.hits = 0
        walk_volume(&vc, self.p1, self.vec, dt_sampler, <void*> &ia)
        if ia.hits > 0:
            return 2 # a box of non-zero volume cannot be inside a ray
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_cell(self, np.float64_t pos[3],
                               np.float64_t dds[3]) nogil:
        # This is terribly inefficient for Octrees.  For grids, it will never
        # get called.
        cdef int i
        cdef np.float64_t left_edge[3]
        cdef np.float64_t right_edge[3]
        for i in range(3):
            left_edge[i] = pos[i] - dds[i]/2.0
            right_edge[i] = pos[i] + dds[i]/2.0
        return self.select_bbox(left_edge, right_edge)

    def _hash_vals(self):
        return (("p1[0]", self.p1[0]),
                ("p1[1]", self.p1[1]),
                ("p1[2]", self.p1[2]),
                ("p2[0]", self.p2[0]),
                ("p2[1]", self.p2[1]),
                ("p2[2]", self.p2[2]),
                ("vec[0]", self.vec[0]),
                ("vec[1]", self.vec[1]),
                ("vec[2]", self.vec[2]))

    def _get_state_attnames(self):
        return ("p1", "p2", "vec")

ray_selector = RaySelector
