cdef class SelectorObject:

    def __cinit__(self, dobj, *args):
        self._hash_initialized = 0
        cdef np.float64_t [:] DLE
        cdef np.float64_t [:] DRE
        min_level = getattr(dobj, "min_level", None)
        max_level = getattr(dobj, "max_level", None)
        if min_level is None:
            min_level = 0
        if max_level is None:
            max_level = 99
        self.min_level = min_level
        self.max_level = max_level
        self.overlap_cells = 0

        ds = getattr(dobj, 'ds', None)
        if ds is None:
            for i in range(3):
                # NOTE that this is not universal.
                self.domain_width[i] = 1.0
                self.periodicity[i] = False
        else:
            DLE = _ensure_code(ds.domain_left_edge)
            DRE = _ensure_code(ds.domain_right_edge)
            for i in range(3):
                self.domain_width[i] = DRE[i] - DLE[i]
                self.domain_center[i] = DLE[i] + 0.5 * self.domain_width[i]
                self.periodicity[i] = ds.periodicity[i]

    def get_periodicity(self):
        cdef int i
        cdef np.ndarray[np.uint8_t, ndim=1] periodicity
        periodicity = np.zeros(3, dtype='uint8')
        for i in range(3):
            periodicity[i] = self.periodicity[i]
        return periodicity

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def select_grids(self,
                     np.ndarray[np.float64_t, ndim=2] left_edges,
                     np.ndarray[np.float64_t, ndim=2] right_edges,
                     np.ndarray[np.int32_t, ndim=2] levels):
        cdef int i, n
        cdef int ng = left_edges.shape[0]
        cdef np.ndarray[np.uint8_t, ndim=1] gridi = np.zeros(ng, dtype='uint8')
        cdef np.float64_t LE[3]
        cdef np.float64_t RE[3]
        _ensure_code(left_edges)
        _ensure_code(right_edges)
        with nogil:
            for n in range(ng):
                # Call our selector function
                # Check if the sphere is inside the grid
                for i in range(3):
                    LE[i] = left_edges[n, i]
                    RE[i] = right_edges[n, i]
                gridi[n] = self.select_grid(LE, RE, levels[n, 0])
        return gridi.astype("bool")

    def count_octs(self, OctreeContainer octree, int domain_id = -1):
        cdef oct_visitors.CountTotalOcts visitor
        visitor = oct_visitors.CountTotalOcts(octree, domain_id)
        octree.visit_all_octs(self, visitor)
        return visitor.index

    def count_oct_cells(self, OctreeContainer octree, int domain_id = -1):
        cdef oct_visitors.CountTotalCells visitor
        visitor = oct_visitors.CountTotalCells(octree, domain_id)
        octree.visit_all_octs(self, visitor)
        return visitor.index

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void recursively_visit_octs(self, Oct *root,
                        np.float64_t pos[3], np.float64_t dds[3],
                        int level,
                        OctVisitor visitor,
                        int visit_covered = 0):
        # visit_covered tells us whether this octree supports partial
        # refinement.  If it does, we need to handle this specially -- first
        # we visit *this* oct, then we make a second pass to check any child
        # octs.
        cdef np.float64_t LE[3]
        cdef np.float64_t RE[3]
        cdef np.float64_t sdds[3]
        cdef np.float64_t spos[3]
        cdef int i, j, k, res
        cdef Oct *ch
        # Remember that pos is the *center* of the oct, and dds is the oct
        # width.  So to get to the edges, we add/subtract half of dds.
        for i in range(3):
            # sdds is the cell width
            sdds[i] = dds[i]/2.0
            LE[i] = pos[i] - dds[i]/2.0
            RE[i] = pos[i] + dds[i]/2.0
        #print(LE[0], RE[0], LE[1], RE[1], LE[2], RE[2])
        res = self.select_grid(LE, RE, level, root)
        if res == 1 and visitor.domain > 0 and root.domain != visitor.domain:
            res = -1
        cdef int increment = 1
        cdef int next_level, this_level
        # next_level: an int that says whether or not we can progress to children
        # this_level: an int that says whether or not we can select from this
        # level
        next_level = this_level = 1
        if res == -1:
            # This happens when we do domain selection but the oct has
            # children.  This would allow an oct to pass to its children but
            # not get accessed itself.
            next_level = 1
            this_level = 0
        elif level == self.max_level:
            next_level = 0
        elif level < self.min_level or level > self.max_level:
            this_level = 0
        if res == 0 and this_level == 1:
            return
        # Now we visit all our children.  We subtract off sdds for the first
        # pass because we center it on the first cell.
        cdef int iter = 1 - visit_covered # 2 if 1, 1 if 0.
        # So the order here goes like so.  If visit_covered is 1, which usually
        # comes from "partial_coverage", we visit the components of a zone even
        # if it has children.  But in general, the first iteration through, we
        # visit each cell.  This means that only if visit_covered is true do we
        # visit potentially covered cells.  The next time through, we visit
        # child cells.
        while iter < 2:
            spos[0] = pos[0] - sdds[0]/2.0
            for i in range(2):
                spos[1] = pos[1] - sdds[1]/2.0
                for j in range(2):
                    spos[2] = pos[2] - sdds[2]/2.0
                    for k in range(2):
                        ch = NULL
                        # We only supply a child if we are actually going to
                        # look at the next level.
                        if root.children != NULL and next_level == 1:
                            ch = root.children[cind(i, j, k)]
                        if iter == 1 and next_level == 1 and ch != NULL:
                            # Note that visitor.pos is always going to be the
                            # position of the Oct -- it is *not* always going
                            # to be the same as the position of the cell under
                            # investigation.
                            visitor.pos[0] = (visitor.pos[0] << 1) + i
                            visitor.pos[1] = (visitor.pos[1] << 1) + j
                            visitor.pos[2] = (visitor.pos[2] << 1) + k
                            visitor.level += 1
                            self.recursively_visit_octs(
                                ch, spos, sdds, level + 1, visitor,
                                visit_covered)
                            visitor.pos[0] = (visitor.pos[0] >> 1)
                            visitor.pos[1] = (visitor.pos[1] >> 1)
                            visitor.pos[2] = (visitor.pos[2] >> 1)
                            visitor.level -= 1
                        elif this_level == 1 and visitor.oref > 0:
                            visitor.global_index += increment
                            increment = 0
                            self.visit_oct_cells(root, ch, spos, sdds,
                                                 visitor, i, j, k)
                        elif this_level == 1 and increment == 1:
                            visitor.global_index += increment
                            increment = 0
                            visitor.ind[0] = visitor.ind[1] = visitor.ind[2] = 0
                            visitor.visit(root, 1)
                        spos[2] += sdds[2]
                    spos[1] += sdds[1]
                spos[0] += sdds[0]
            this_level = 0 # We turn this off for the second pass.
            iter += 1

    cdef void visit_oct_cells(self, Oct *root, Oct *ch,
                              np.float64_t spos[3], np.float64_t sdds[3],
                              OctVisitor visitor, int i, int j, int k):
        # We can short-circuit the whole process if data.oref == 1.
        # This saves us some funny-business.
        cdef int selected
        if visitor.oref == 1:
            selected = self.select_cell(spos, sdds)
            if ch != NULL:
                selected *= self.overlap_cells
            # visitor.ind refers to the cell, not to the oct.
            visitor.ind[0] = i
            visitor.ind[1] = j
            visitor.ind[2] = k
            visitor.visit(root, selected)
            return
        # Okay, now that we've got that out of the way, we have to do some
        # other checks here.  In this case, spos[] is the position of the
        # center of a *possible* oct child, which means it is the center of a
        # cluster of cells.  That cluster might have 1, 8, 64, ... cells in it.
        # But, we can figure it out by calculating the cell dds.
        cdef np.float64_t dds[3]
        cdef np.float64_t pos[3]
        cdef int ci, cj, ck
        cdef int nr = (1 << (visitor.oref - 1))
        for ci in range(3):
            dds[ci] = sdds[ci] / nr
        # Boot strap at the first index.
        pos[0] = (spos[0] - sdds[0]/2.0) + dds[0] * 0.5
        for ci in range(nr):
            pos[1] = (spos[1] - sdds[1]/2.0) + dds[1] * 0.5
            for cj in range(nr):
                pos[2] = (spos[2] - sdds[2]/2.0) + dds[2] * 0.5
                for ck in range(nr):
                    selected = self.select_cell(pos, dds)
                    if ch != NULL:
                        selected *= self.overlap_cells
                    visitor.ind[0] = ci + i * nr
                    visitor.ind[1] = cj + j * nr
                    visitor.ind[2] = ck + k * nr
                    visitor.visit(root, selected)
                    pos[2] += dds[2]
                pos[1] += dds[1]
            pos[0] += dds[0]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_grid(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3],
                               np.int32_t level, Oct *o = NULL) nogil:
        if level < self.min_level or level > self.max_level: return 0
        return self.select_bbox(left_edge, right_edge)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_grid_edge(self, np.float64_t left_edge[3],
                                    np.float64_t right_edge[3],
                                    np.int32_t level, Oct *o = NULL) nogil:
        if level < self.min_level or level > self.max_level: return 0
        return self.select_bbox_edge(left_edge, right_edge)

    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        return 0

    cdef int select_point(self, np.float64_t pos[3]) nogil:
        return 0

    cdef int select_sphere(self, np.float64_t pos[3], np.float64_t radius) nogil:
        return 0

    cdef int select_bbox(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        """
        Returns:
          0: If the selector does not touch the bounding box.
          1: If the selector overlaps the bounding box anywhere.
        """
        return 0

    cdef int select_bbox_edge(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        """
        Returns:
          0: If the selector does not touch the bounding box.
          1: If the selector contains the entire bounding box.
          2: If the selector contains part of the bounding box.
        """
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef np.float64_t periodic_difference(self, np.float64_t x1, np.float64_t x2, int d) nogil:
        # domain_width is already in code units, and we assume what is fed in
        # is too.
        cdef np.float64_t rel = x1 - x2
        if self.periodicity[d]:
            if rel > self.domain_width[d] * 0.5:
                rel -= self.domain_width[d]
            elif rel < -self.domain_width[d] * 0.5:
                rel += self.domain_width[d]
        return rel

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def fill_mesh_mask(self, mesh):
        cdef np.float64_t pos[3]
        cdef np.ndarray[np.int64_t, ndim=2] indices
        cdef np.ndarray[np.float64_t, ndim=2] coords
        cdef np.ndarray[np.uint8_t, ndim=1] mask
        cdef int i, j, k, selected
        cdef int npoints, nv = mesh._connectivity_length
        cdef int total = 0
        cdef int offset = mesh._index_offset
        coords = _ensure_code(mesh.connectivity_coords)
        indices = mesh.connectivity_indices
        npoints = indices.shape[0]
        mask = np.zeros(npoints, dtype='uint8')
        for i in range(npoints):
            selected = 0
            for j in range(nv):
                for k in range(3):
                    pos[k] = coords[indices[i, j] - offset, k]
                selected = self.select_point(pos)
                if selected == 1: break
            total += selected
            mask[i] = selected
        if total == 0: return None
        return mask.astype("bool")

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def fill_mesh_cell_mask(self, mesh):
        cdef np.float64_t pos
        cdef np.float64_t le[3]
        cdef np.float64_t re[3]
        cdef np.ndarray[np.int64_t, ndim=2] indices
        cdef np.ndarray[np.float64_t, ndim=2] coords
        cdef np.ndarray[np.uint8_t, ndim=1] mask
        cdef int i, j, k, selected
        cdef int npoints, nv = mesh._connectivity_length
        cdef int ndim = mesh.connectivity_coords.shape[1]
        cdef int total = 0
        cdef int offset = mesh._index_offset
        coords = _ensure_code(mesh.connectivity_coords)
        indices = mesh.connectivity_indices
        npoints = indices.shape[0]
        mask = np.zeros(npoints, dtype='uint8')
        for i in range(npoints):
            selected = 0
            for k in range(3):
                le[k] = 1e60
                re[k] = -1e60
            for j in range(nv):
                for k in range(ndim):
                    pos = coords[indices[i, j] - offset, k]
                    le[k] = fmin(pos, le[k])
                    re[k] = fmax(pos, re[k])
                for k in range(2, ndim - 1, -1):
                    le[k] = self.domain_center[k]
                    re[k] = self.domain_center[k]
            selected = self.select_bbox(le, re)
            total += selected
            mask[i] = selected
        if total == 0: return None
        return mask.astype("bool")

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def fill_mask(self, gobj):
        cdef np.ndarray[np.uint8_t, ndim=3, cast=True] child_mask
        child_mask = gobj.child_mask
        cdef np.ndarray[np.uint8_t, ndim=3] mask
        cdef int dim[3]
        _ensure_code(gobj.dds)
        _ensure_code(gobj.LeftEdge)
        _ensure_code(gobj.RightEdge)
        cdef np.ndarray[np.float64_t, ndim=1] odds = gobj.dds.d
        cdef np.ndarray[np.float64_t, ndim=1] oleft_edge = gobj.LeftEdge.d
        cdef np.ndarray[np.float64_t, ndim=1] oright_edge = gobj.RightEdge.d
        cdef int i
        cdef np.float64_t dds[3]
        cdef np.float64_t left_edge[3]
        cdef np.float64_t right_edge[3]
        for i in range(3):
            dds[i] = odds[i]
            dim[i] = gobj.ActiveDimensions[i]
            left_edge[i] = oleft_edge[i]
            right_edge[i] = oright_edge[i]
        mask = np.zeros(gobj.ActiveDimensions, dtype='uint8')
        # Check for the level bounds
        cdef np.int32_t level = gobj.Level
        # We set this to 1 if we ignore child_mask
        cdef int total
        total = self.fill_mask_selector(left_edge, right_edge, dds, dim,
                                        child_mask, mask, level)
        if total == 0: return None
        return mask.astype("bool")

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int fill_mask_selector(self, np.float64_t left_edge[3],
                                np.float64_t right_edge[3],
                                np.float64_t dds[3], int dim[3],
                                np.ndarray[np.uint8_t, ndim=3, cast=True] child_mask,
                                np.ndarray[np.uint8_t, ndim=3] mask,
                                int level):
        cdef int i, j, k
        cdef int total = 0, this_level = 0
        cdef np.float64_t pos[3]
        if level < self.min_level or level > self.max_level:
            return 0
        if level == self.max_level:
            this_level = 1
        with nogil:
            for i in range(dim[0]):
                pos[0] = left_edge[0] + (i + 0.5) * dds[0]
                for j in range(dim[1]):
                    pos[1] = left_edge[1] + (j + 0.5) * dds[1]
                    for k in range(dim[2]):
                        pos[2] = left_edge[2] + (k + 0.5) * dds[2]
                        if child_mask[i, j, k] == 1 or this_level == 1:
                            mask[i, j, k] = self.select_cell(pos, dds)
                            total += mask[i, j, k]
        return total

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void visit_grid_cells(self, GridVisitorData *data,
                              grid_visitor_function *func,
                              np.uint8_t *cached_mask = NULL):
        # This function accepts a grid visitor function, the data that
        # corresponds to the current grid being examined (the most important
        # aspect of which is the .grid attribute, along with index values and
        # void* pointers to arrays) and a possibly-pre-generated cached mask.
        # Each cell is visited with the grid visitor function.
        cdef np.float64_t left_edge[3]
        cdef np.float64_t right_edge[3]
        cdef np.float64_t dds[3]
        cdef int dim[3]
        cdef int this_level = 0, level, i
        cdef np.float64_t pos[3]
        level = data.grid.level
        if level < self.min_level or level > self.max_level:
            return
        if level == self.max_level:
            this_level = 1
        cdef np.uint8_t child_masked, selected
        for i in range(3):
            left_edge[i] = data.grid.left_edge[i]
            right_edge[i] = data.grid.right_edge[i]
            dds[i] = (right_edge[i] - left_edge[i])/data.grid.dims[i]
            dim[i] = data.grid.dims[i]
        with nogil:
            pos[0] = left_edge[0] + dds[0] * 0.5
            data.pos[0] = 0
            for i in range(dim[0]):
                pos[1] = left_edge[1] + dds[1] * 0.5
                data.pos[1] = 0
                for j in range(dim[1]):
                    pos[2] = left_edge[2] + dds[2] * 0.5
                    data.pos[2] = 0
                    for k in range(dim[2]):
                        # We short-circuit if we have a cache; if we don't, we
                        # only set selected to true if it's *not* masked by a
                        # child and it *is* selected.
                        if cached_mask != NULL:
                            selected = ba_get_value(cached_mask,
                                                    data.global_index)
                        else:
                            if this_level == 1:
                                child_masked = 0
                            else:
                                child_masked = check_child_masked(data)
                            if child_masked == 0:
                                selected = self.select_cell(pos, dds)
                            else:
                                selected = 0
                        func(data, selected)
                        data.global_index += 1
                        pos[2] += dds[2]
                        data.pos[2] += 1
                    pos[1] += dds[1]
                    data.pos[1] += 1
                pos[0] += dds[0]
                data.pos[0] += 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def count_points(self, np.ndarray[floating, ndim=1] x,
                           np.ndarray[floating, ndim=1] y,
                           np.ndarray[floating, ndim=1] z,
                           radii):
        cdef int count = 0
        cdef int i
        cdef np.float64_t pos[3]
        cdef np.float64_t radius
        cdef np.float64_t[:] _radii
        if radii is not None:
            _radii = np.atleast_1d(np.array(radii, dtype='float64'))
        else:
            _radii = np.array([0.0], dtype='float64')
        _ensure_code(x)
        _ensure_code(y)
        _ensure_code(z)
        with nogil:
            for i in range(x.shape[0]):
                pos[0] = x[i]
                pos[1] = y[i]
                pos[2] = z[i]
                if _radii.shape[0] == 1:
                    radius = _radii[0]
                else:
                    radius = _radii[i]
                if radius == 0:
                    count += self.select_point(pos)
                else:
                    count += self.select_sphere(pos, radius)
        return count

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def select_points(self,
                      np.ndarray[floating, ndim=1] x,
                      np.ndarray[floating, ndim=1] y,
                      np.ndarray[floating, ndim=1] z,
                      radii):
        cdef int count = 0
        cdef int i
        cdef np.float64_t pos[3]
        cdef np.float64_t radius
        cdef np.ndarray[np.uint8_t, ndim=1] mask
        cdef np.float64_t[:] _radii
        if radii is not None:
            _radii = np.atleast_1d(np.array(radii, dtype='float64'))
        else:
            _radii = np.array([0.0], dtype='float64')
        mask = np.empty(x.shape[0], dtype='uint8')
        _ensure_code(x)
        _ensure_code(y)
        _ensure_code(z)


        # this is to allow selectors to optimize the point vs
        # 0-radius sphere case.  These two may have different
        # effects for 0-volume selectors, however (collision
        # between a ray and a point is null, while ray and a
        # sphere is allowed)
        with nogil:
            for i in range(x.shape[0]) :
                pos[0] = x[i]
                pos[1] = y[i]
                pos[2] = z[i]
                if _radii.shape[0] == 1:
                    radius = 0
                else:
                    radius = _radii[i]
                if radius == 0:
                    mask[i] = self.select_point(pos)
                else:
                    mask[i] = self.select_sphere(pos, radius)
                count += mask[i]
        if count == 0: return None
        return mask.view("bool")

    def __hash__(self):
        # convert data to be hashed to a byte array, which FNV algorithm expects
        if self._hash_initialized == 1:
            return self._hash
        hash_data = bytearray()
        for v in self._hash_vals() + self._base_hash():
            if isinstance(v, tuple):
                hash_data.extend(v[0].encode('ascii'))
                hash_data.extend(repr(v[1]).encode('ascii'))
            else:
                hash_data.extend(repr(v).encode('ascii'))
        cdef np.int64_t hash_value = fnv_hash(hash_data)
        self._hash = hash_value
        self._hash_initialized = 1
        return hash_value

    def _hash_vals(self):
        raise NotImplementedError

    def _base_hash(self):
        return (("min_level", self.min_level),
                ("max_level", self.max_level),
                ("overlap_cells", self.overlap_cells),
                ("periodicity[0]", self.periodicity[0]),
                ("periodicity[1]", self.periodicity[1]),
                ("periodicity[2]", self.periodicity[2]),
                ("domain_width[0]", self.domain_width[0]),
                ("domain_width[1]", self.domain_width[1]),
                ("domain_width[2]", self.domain_width[2]))

    def _get_state_attnames(self):
        # return a tupe of attr names for __setstate__: implement for each subclass
        raise NotImplementedError

    def __getstate__(self):
        # returns a tuple containing (attribute name, attribute value) tuples needed to
        # rebuild the state:
        base_atts = ("min_level", "max_level", "overlap_cells",
                     "periodicity", "domain_width", "domain_center")
        child_atts = self._get_state_attnames()

        # assemble the state_tuple (('a1', a1val), ('a2', a2val),...)
        state_tuple = ()
        for fld in base_atts + child_atts:
            state_tuple += ((fld, getattr(self, fld)), )
        return state_tuple

    def __getnewargs__(self):
        # __setstate__ will always call __cinit__, this pickle hook returns arguments
        # to __cinit__. We will give it None so we dont error then set attributes in
        # __setstate__ Note that we could avoid this by making dobj an optional argument
        # to __cinit__
        return (None, )

    def __setstate__(self, state_tuple):
        # parse and set attributes from the state_tuple: (('a1',a1val),('a2',a2val),...)
        for attr in state_tuple:
            setattr(self, attr[0], attr[1])
