cdef class CutRegionSelector(SelectorObject):
    cdef set _positions
    cdef tuple _conditionals

    def __init__(self, dobj):
        axis_name = dobj.ds.coordinates.axis_name
        positions = np.array([dobj['index', axis_name[0]],
                              dobj['index', axis_name[1]],
                              dobj['index', axis_name[2]]]).T
        self._conditionals = tuple(dobj.conditionals)
        self._positions = set(tuple(position) for position in positions)

    cdef int select_bbox(self,  np.float64_t left_edge[3],
                     np.float64_t right_edge[3]) nogil:
        return 1

    cdef int select_bbox_dge(self,  np.float64_t left_edge[3],
                     np.float64_t right_edge[3]) nogil:
        return 1

    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        with gil:
            if (pos[0], pos[1], pos[2]) in self._positions:
                return 1
            else:
                return 0

    cdef int select_point(self, np.float64_t pos[3]) nogil:
        return 1

    cdef int select_sphere(self, np.float64_t pos[3], np.float64_t radius) nogil:
        return 1

    def _hash_vals(self):
        t = ()
        for i, c in enumerate(self._conditionals):
            t += ("conditional[%s]" % i, c)
        return ("conditionals", t)

cut_region_selector = CutRegionSelector
