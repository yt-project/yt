"""
Geometry selection routine imports.




"""


cimport numpy as np
from grid_visitors cimport (
    GridTreeNode,
    GridVisitorData,
    check_child_masked,
    grid_visitor_function,
)
from oct_container cimport OctreeContainer
from oct_visitors cimport Oct, OctVisitor

from yt.utilities.lib.fp_utils cimport _ensure_code
from yt.utilities.lib.geometry_utils cimport decode_morton_64bit


cdef class SelectorObject:
    cdef public np.int32_t min_level
    cdef public np.int32_t max_level
    cdef public int overlap_cells
    cdef public np.float64_t domain_width[3]
    cdef public np.float64_t domain_center[3]
    cdef public bint periodicity[3]
    cdef bint _hash_initialized
    cdef np.int64_t _hash

    cdef void recursively_visit_octs(self, Oct *root,
                        np.float64_t pos[3], np.float64_t dds[3],
                        int level,
                        OctVisitor visitor,
                        int visit_covered = ?)
    cdef void visit_oct_cells(self, Oct *root, Oct *ch,
                              np.float64_t spos[3], np.float64_t sdds[3],
                              OctVisitor visitor, int i, int j, int k)
    cdef int select_grid(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3],
                               np.int32_t level, Oct *o = ?) nogil
    cdef int select_grid_edge(self, np.float64_t left_edge[3],
                                    np.float64_t right_edge[3],
                                    np.int32_t level, Oct *o = ?) nogil
    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil

    cdef int select_point(self, np.float64_t pos[3]) nogil
    cdef int select_sphere(self, np.float64_t pos[3], np.float64_t radius) nogil
    cdef int select_bbox(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil
    cdef int select_bbox_edge(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil
    cdef int fill_mask_selector_regular_grid(self, np.float64_t left_edge[3],
                                             np.float64_t right_edge[3],
                                             np.float64_t dds[3], int dim[3],
                                             np.ndarray[np.uint8_t, ndim=3, cast=True] child_mask,
                                             np.ndarray[np.uint8_t, ndim=3] mask,
                                             int level)
    cdef int fill_mask_selector(self, np.float64_t left_edge[3],
                                np.float64_t right_edge[3],
                                np.float64_t **dds, int dim[3],
                                np.ndarray[np.uint8_t, ndim=3, cast=True] child_mask,
                                np.ndarray[np.uint8_t, ndim=3] mask,
                                int level)
    cdef void visit_grid_cells(self, GridVisitorData *data,
                    grid_visitor_function *func, np.uint8_t *cached_mask = ?)

    # compute periodic distance (if periodicity set)
    # assuming 0->domain_width[d] coordinates
    cdef np.float64_t periodic_difference(
        self, np.float64_t x1, np.float64_t x2, int d) nogil

cdef class AlwaysSelector(SelectorObject):
    pass

cdef class OctreeSubsetSelector(SelectorObject):
    cdef public SelectorObject base_selector
    cdef public np.int64_t domain_id

cdef class BooleanSelector(SelectorObject):
    cdef public SelectorObject sel1
    cdef public SelectorObject sel2

cdef inline np.float64_t _periodic_dist(np.float64_t x1, np.float64_t x2,
                                        np.float64_t dw, bint periodic) nogil:
    cdef np.float64_t rel = x1 - x2
    if not periodic: return rel
    if rel > dw * 0.5:
        rel -= dw
    elif rel < -dw * 0.5:
        rel += dw
    return rel
