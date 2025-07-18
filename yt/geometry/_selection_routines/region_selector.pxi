cdef class RegionSelector(SelectorObject):
    cdef public np.float64_t left_edge[3]
    cdef public np.float64_t right_edge[3]
    cdef public np.float64_t right_edge_shift[3]
    cdef public bint is_all_data
    cdef public bint loose_selection
    cdef public bint check_period[3]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def __init__(self, dobj):
        cdef int i
        # We are modifying dobj.left_edge and dobj.right_edge , so here we will
        # do an in-place conversion of those arrays.
        cdef np.float64_t[:] RE = dobj.right_edge.copy()
        cdef np.float64_t[:] LE = dobj.left_edge.copy()
        cdef const np.float64_t[:] DW = dobj.ds.domain_width
        cdef const np.float64_t[:] DLE = dobj.ds.domain_left_edge
        cdef const np.float64_t[:] DRE = dobj.ds.domain_right_edge
        le_all = (np.array(LE) == dobj.ds.domain_left_edge).all()
        re_all = (np.array(RE) == dobj.ds.domain_right_edge).all()
        # If we have a bounding box, then we should *not* revert to all data
        domain_override = getattr(dobj.ds, '_domain_override', False)
        if le_all and re_all and not domain_override:
            self.is_all_data = True
        else:
            self.is_all_data = False
        cdef np.float64_t region_width[3]
        cdef bint p[3]
        # This is for if we want to include zones that overlap and whose
        # centers are not strictly included.
        self.loose_selection = getattr(dobj, "loose_selection", False)

        for i in range(3):
            self.check_period[i] = False
            region_width[i] = RE[i] - LE[i]
            p[i] = dobj.ds.periodicity[i]
            if region_width[i] <= 0:
                raise RuntimeError(
                    "Region right edge[%s] < left edge: width = %s" % (
                        i, region_width[i]))

        for i in range(3):

            if p[i]:
                # First, we check if any criteria requires a period check,
                # without any adjustments.  This is for short-circuiting the
                # short-circuit of the loop down below in mask filling.
                if LE[i] < DLE[i] or LE[i] > DRE[i] or RE[i] > DRE[i]:
                    self.check_period[i] = True
                # shift so left_edge guaranteed in domain
                if LE[i] < DLE[i]:
                    LE[i] += DW[i]
                    RE[i] += DW[i]
                elif LE[i] > DRE[i]:
                    LE[i] -= DW[i]
                    RE[i] -= DW[i]
            else:
                if LE[i] < DLE[i] or RE[i] > DRE[i]:
                    raise RuntimeError(
                        "yt attempted to read outside the boundaries of "
                        "a non-periodic domain along dimension %s.\n"
                        "Region left edge = %s, Region right edge = %s\n"
                        "Dataset left edge = %s, Dataset right edge = %s\n\n"
                        "This commonly happens when trying to compute ghost cells "
                        "up to the domain boundary. Two possible solutions are to "
                        "select a smaller region that does not border domain edge "
                        "(see https://yt-project.org/docs/analyzing/objects.html?highlight=region)\n"
                        "or override the periodicity with\n"
                        "ds.force_periodicity()" % \
                        (i, dobj.left_edge[i], dobj.right_edge[i],
                         dobj.ds.domain_left_edge[i], dobj.ds.domain_right_edge[i])
                    )
            # Already ensured in code
            self.left_edge[i] = LE[i]
            self.right_edge[i] = RE[i]
            self.right_edge_shift[i] = RE[i] - DW[i]
            if not self.periodicity[i]:
                self.right_edge_shift[i] = -np.inf

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_bbox(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) noexcept nogil:
        cdef int i
        for i in range(3):
            if (right_edge[i] < self.left_edge[i] and \
                left_edge[i] >= self.right_edge_shift[i]) or \
                left_edge[i] >= self.right_edge[i]:
                return 0
        return 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_bbox_edge(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) noexcept nogil:
        cdef int i
        for i in range(3):
            if (right_edge[i] < self.left_edge[i] and \
                left_edge[i] >= self.right_edge_shift[i]) or \
                left_edge[i] >= self.right_edge[i]:
                return 0
        for i in range(3):
            if left_edge[i] < self.right_edge_shift[i]:
                if right_edge[i] >= self.right_edge_shift[i]:
                    return 2
            elif left_edge[i] < self.left_edge[i] or \
                 right_edge[i] >= self.right_edge[i]:
                return 2
        return 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) noexcept nogil:
        cdef np.float64_t left_edge[3]
        cdef np.float64_t right_edge[3]
        cdef int i
        if self.loose_selection:
            for i in range(3):
                left_edge[i] = pos[i] - dds[i]*0.5
                right_edge[i] = pos[i] + dds[i]*0.5
            return self.select_bbox(left_edge, right_edge)
        return self.select_point(pos)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_point(self, np.float64_t pos[3]) noexcept nogil:
        cdef int i
        for i in range(3):
            if (self.right_edge_shift[i] <= pos[i] < self.left_edge[i]) or \
               pos[i] >= self.right_edge[i]:
                return 0
        return 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_sphere(self, np.float64_t pos[3], np.float64_t radius) noexcept nogil:
        # adapted from http://stackoverflow.com/a/4579192/1382869
        cdef int i
        cdef np.float64_t p
        cdef np.float64_t r2 = radius**2
        cdef np.float64_t dmin = 0
        cdef np.float64_t d = 0
        for i in range(3):
            if (pos[i]+radius < self.left_edge[i] and \
                pos[i]-radius >= self.right_edge_shift[i]):
                d = self.periodic_difference(pos[i], self.left_edge[i], i)
            elif pos[i]-radius > self.right_edge[i]:
                d = self.periodic_difference(pos[i], self.right_edge[i], i)
            dmin += d*d
        return int(dmin <= r2)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int fill_mask_selector_regular_grid(self, np.float64_t left_edge[3],
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
        cdef np.int64_t si[3]
        cdef np.int64_t ei[3]
        for i in range(3):
            if not self.check_period[i]:
                si[i] = <np.int64_t> ((self.left_edge[i] - left_edge[i])/dds[i])
                ei[i] = <np.int64_t> ((self.right_edge[i] - left_edge[i])/dds[i])
                si[i] = iclip(si[i] - 1, 0, dim[i])
                ei[i] = iclip(ei[i] + 1, 0, dim[i])
            else:
                si[i] = 0
                ei[i] = dim[i]
        with nogil:

            for i in range(si[0], ei[0]):
                pos[0] = left_edge[0] + (i + 0.5) * dds[0]
                for j in range(si[1], ei[1]):
                    pos[1] = left_edge[1] + (j + 0.5) * dds[1]
                    for k in range(si[2], ei[2]):
                        pos[2] = left_edge[2] + (k + 0.5) * dds[2]
                        if child_mask[i, j, k] == 1 or this_level == 1:
                            mask[i, j, k] = self.select_cell(pos, dds)
                            total += mask[i, j, k]
        return total

    def _hash_vals(self):
        return (("left_edge[0]", self.left_edge[0]),
                ("left_edge[1]", self.left_edge[1]),
                ("left_edge[2]", self.left_edge[2]),
                ("right_edge[0]", self.right_edge[0]),
                ("right_edge[1]", self.right_edge[1]),
                ("right_edge[2]", self.right_edge[2]))

    def _get_state_attnames(self):
        return ('left_edge', 'right_edge', 'right_edge_shift', 'check_period',
                'is_all_data', 'loose_selection')

region_selector = RegionSelector
