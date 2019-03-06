cimport numpy as np
cimport cython

from yt.utilities.lib.cykdtree.kdtree cimport KDTree, uint64_t
from yt.utilities.lib.bounded_priority_queue cimport BoundedPriorityQueue

cdef struct axes_range:
    int start
    int stop
    int step

cdef int set_axes_range(axes_range *axes, int skipaxis)

cdef int find_neighbors(np.float64_t * pos, np.float64_t[:, ::1] tree_positions,
                        BoundedPriorityQueue queue, KDTree * c_tree,
                        uint64_t skipidx, axes_range * axes) nogil except -1
